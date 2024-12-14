#include <cmath>
#include <ctime>
#include <highfive/H5File.hpp>
#include <initializer_list>
#include <limits>
#include <list>
#include <random>
#include <map>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <nanoflann.hpp>

#include <GeographicLib/Geodesic.hpp>

#include <gswteos-10.h>

using namespace std;
using namespace GeographicLib;
using namespace nlohmann;

enum class Cols: int {
	depth=0,
	latitude=1,
	longitude=2,
	so=3,
	thetao=4,
	uo=5,
	vo=6,
	num_column=7
};

struct Cell {
	bool exists=false;
	double density,u,v;
};

constexpr double eps = 1e-7;
struct EpsCmp {
	bool operator()(double x, double y) const {
		return abs(x-y)>eps && x<y;
	}
};

struct Parameters {
	double lat,lng,depth,east,north,down,mass,volume;
	double drag, area;
};

double sgn(double x) {
	return x>eps ? 1 : x < -eps ? -1 : 0;
}

constexpr double inf = numeric_limits<double>::infinity();

struct BBox {
	array<double,3> b_min,b_max;
	BBox() {
		b_min.fill(inf), b_max.fill(-inf);
	}
	void add(array<double,3> pos) {
		for (int i=0; i<3; i++) b_min[i]=min(pos[i], b_min[i]);
		for (int i=0; i<3; i++) b_max[i]=max(pos[i], b_max[i]);
	}
};

struct Data {
	vector<double> depths, lats, lngs;
	map<double, int, EpsCmp> lat_idx, lng_idx, depth_idx;
	vector<Cell> cells;
	int num_lat, num_lng, num_depth;

	Geodesic const& geo = Geodesic::WGS84();

	Data(string const& path) {
		HighFive::File file(path);
		auto dataset = file.getDataSet("data");
		auto data = dataset.read<vector<array<double,unsigned(Cols::num_column)>>>();

		for (auto& v: data) {
			depth_idx.insert({v[int(Cols::depth)],-1});
			lat_idx.insert({v[int(Cols::latitude)],-1});
			lng_idx.insert({v[int(Cols::longitude)],-1});
		}

		for (auto [a,_]: depth_idx) depths.push_back(a);
		for (auto [a,_]: lat_idx) lats.push_back(a);
		for (auto [a,_]: lng_idx) lngs.push_back(a);

		for (auto map: {&lat_idx, &lng_idx, &depth_idx}) {
			int i=0;
			for (auto& [a,b]: *map) b=i++;

			map->insert({-inf, -1});
			map->insert({inf, -1});
		}

		num_lat=lats.size(), num_lng=lngs.size(), num_depth=depths.size();
		cells.resize(num_lat*num_lng*num_depth);

		for (auto& v: data) {
			bool nan=false;
			for (double x: v) {
				if (isnan(x)) {nan=true; break;}
			}
			if (nan) continue;

			double depth = v[int(Cols::depth)];
			double lat = v[int(Cols::latitude)];
			double lng = v[int(Cols::longitude)];

			double practical_salinity = v[int(Cols::so)];
			double potential_temp = v[int(Cols::thetao)];

			// assume both dynamic height anomaly and sea surface geopotential are zero
			// i don't really know what those are but they don't sound like things that greatly affect the pressure at this depth...
			double pressure = gsw_p_from_z(-depth, lat, 0,0);
			double absolute_salinity = gsw_sa_from_sp(practical_salinity, pressure, lng, lat);
			double conservative_temperature = gsw_ct_from_pt(absolute_salinity, potential_temp);
			double density = gsw_rho(absolute_salinity, conservative_temperature, pressure);
			if (isnan(density)) {
				throw runtime_error("nan density!");
			}

			at(lat_idx[lat],lng_idx[lng],depth_idx[depth]) = Cell {
				.exists=true, .density=density,
				.u=v[int(Cols::uo)], .v=v[int(Cols::vo)]
			};
		}
	}

	Cell& at(int lat_i, int lng_i, int depth_i) {
		return cells[lat_i + num_lat*(lng_i + num_lng*depth_i)];
	};

	bool at_interp(int depth_i, int lat_i, int lng_i, double cur_lat, double cur_lng, double cur_depth, Cell& interpolated) {
		double tot = (depths[depth_i]-depths[depth_i-1])
			*(lats[lat_i]-lats[lat_i-1])
			*(lngs[lng_i]-lngs[lng_i-1]);

		for (int z=0; z<=1; z++) {
			for (int y=0; y<=1; y++) {
				for (int x=0; x<=1; x++) {
					auto& cell = at(lat_i-1+x, lng_i-1+y, depth_i-1+z);
					if (!cell.exists) return true;

					double coeff = (z ? depths[depth_i]-cur_depth : cur_depth-depths[depth_i-1])
						*(x ? lats[lat_i]-cur_lat : cur_lat-lats[lat_i-1])
						*(y ? lngs[lng_i]-cur_lng : cur_lng-lngs[lng_i-1])/tot;

					interpolated.density += coeff*cell.density;
					interpolated.u += coeff*cell.u;
					interpolated.v += coeff*cell.v;
				}
			}
		}

		return false;
	};

	Cell get_cell(double lat, double lng, double depth) {
		int depth_i = depth_idx.upper_bound(depth)->second;
		int lat_i = lat_idx.upper_bound(lat)->second;
		int lng_i = lng_idx.upper_bound(lng)->second;

		for (int x: {depth_i, lat_i, lng_i}) {
			if (x==-1) return {.exists=false};
		}

		Cell out {.exists=true};
		if (at_interp(depth_i,lat_i,lng_i,lat,lng,depth,out))
			return {.exists=false};
		return out;
	}
};

struct Point {
	double lat,lng,depth,north,east,down;
	Point operator+(Point const& o) const {
		return Point {
			lat+o.lat,lng+o.lng,depth+o.depth,
			north+o.north,east+o.east,down+o.down
		};
	}
	Point operator-(Point const& o) const {
		return Point {
			lat-o.lat,lng-o.lng,depth-o.depth,
			north-o.north,east-o.east,down-o.down
		};
	}
	Point operator/(double o) const {
		return Point {
			lat/o,lng/o,depth/o,
			north/o,east/o,down/o
		};
	}
	Point operator*(double o) const {
		return Point {
			lat*o,lng*o,depth*o,
			north*o,east*o,down*o
		};
	}
};

struct Simulation {
	vector<vector<Point>> points;
};

constexpr double dt = 5;
constexpr double maxt = 24*3600;
constexpr int log_every = 20;
constexpr double log_every_t = log_every*dt;
constexpr int num_sim=1000;

Simulation simulate(Data& data, Parameters params) {
	constexpr double noise_update_second = 0.001;
	constexpr double g = 9.81;

	Simulation s;
	ofstream out("./out.csv");

	minstd_rand rng(random_device{}());
	normal_distribution<> u_dist(0,0.1);
	normal_distribution<> down_dist(0,0.01);
	normal_distribution<> v_dist(0,0.1);

	out<<"t,lat,lng,depth,east,north\n"<<fixed<<setprecision(6);

	auto iter = [](int& idx, vector<double> const& vs, double target) {
		while (idx<vs.size()) {
			if (vs[idx]<target) idx++;
			else if (idx>0 && vs[idx-1]>target) idx--;
			else break;
		}
		
		return idx<=0 || idx>=vs.size();
	};

	BBox box;

	for (int sim_i=0; sim_i<num_sim; sim_i++) {
		auto& sim_points = s.points.emplace_back();

		double t=0;
		int i=0;
		bool crashed=false, exited=false;

		double cur_lat = params.lat, cur_lng = params.lng, cur_depth = params.depth;
		double cur_down = params.down, cur_east = params.east, cur_north = params.north;

		double down_noise=down_dist(rng), u_noise=u_dist(rng), v_noise=v_dist(rng);

		int depth_i = data.depth_idx.upper_bound(cur_depth)->second;
		int lat_i = data.lat_idx.upper_bound(cur_lat)->second;
		int lng_i = data.lng_idx.upper_bound(cur_lng)->second;

		while (t<maxt) {
			if ((i++)%log_every == 0) {
				initializer_list<double> props = {t,cur_lat,cur_lng,cur_depth,cur_east,cur_north};
				for (int prop_i=0; prop_i<props.size(); prop_i++) {
					out<<*(props.begin()+prop_i)<<",\n"[prop_i==props.size()-1];
				}

				sim_points.push_back(Point {
					cur_lat,cur_lng,cur_depth,
					cur_north,cur_east,cur_down
				});

				box.add({cur_lat,cur_lng,cur_depth});
			}

			t+=dt;

			cur_depth += dt*cur_down;

			double ang = atan2(cur_east, cur_north)*180/M_PI;
			if (!isnan(ang)) {
				if (ang<0) ang=360-ang;
				data.geo.Direct(cur_lat, cur_lng, ang, dt*hypot(cur_east, cur_north), cur_lat, cur_lng);
			}

			if (iter(depth_i, data.depths, cur_depth)) {
				crashed=true; break;
			}
			if (iter(lat_i, data.lats, cur_lat)
				|| iter(lng_i, data.lngs, cur_lng)) {
				exited=true; break;
			}

			double a = noise_update_second*dt, b=1-a;
			down_noise = a*down_dist(rng) + b*down_noise;
			u_noise = a*u_dist(rng) + b*u_noise;
			v_noise = a*v_dist(rng) + b*v_noise;

			Cell interpolated = {.density=0, .u=u_noise, .v=v_noise};
			if (data.at_interp(depth_i,lat_i,lng_i,cur_lat,cur_lng,cur_depth,interpolated)) {
				crashed=true;
				break;
			}

			double mul = dt/params.mass;
			double drag = hypot(down_noise - cur_down, interpolated.v - cur_north, interpolated.u - cur_east);

			cur_down += mul*(
				g*(params.mass - interpolated.density*params.volume)
				+ interpolated.density*(down_noise - cur_down)*drag*params.drag*params.area/2
			);

			cur_north += mul*(
				interpolated.density*(interpolated.v - cur_north)*drag*params.drag*params.area/2
			);

			cur_east += mul*(
				interpolated.density*(interpolated.u - cur_east)*drag*params.drag*params.area/2
			);
		}

		out<<"\n";

		cout<<"done simulating #"<<sim_i<<" to t="<<t<<"\n";
		if (crashed) cout<<"crashed into ocean floor or ran aground\n";
		else if (exited) cout<<"exited bounds\n";
	}

	cout<<"outputting cells in bounding box\n";

	ofstream cell_out("./cells.csv");
	cell_out<<"lat,lng,depth,north,east,density\n";

	int depth_i=0, lat_i=0, lng_i=0;
	constexpr double depth_inc = 100, lat_inc=0.01, lng_inc=lat_inc/3;

	for (double depth=max(box.b_min[2], data.depths[0]); depth<=min(box.b_max[2], data.depths.back()); depth+=depth_inc) {
		iter(depth_i, data.depths, depth);
		for (double lat=max(box.b_min[0], data.lats[0]); lat<=min(box.b_max[0], data.lats.back()); lat+=lat_inc) {
			iter(lat_i, data.lats, lat);
			for (double lng=max(box.b_min[1], data.lngs[0]); lng<=min(box.b_max[1], data.lngs.back()); lng+=lng_inc) {
				iter(lng_i, data.lngs, lng);
				Cell interp {0,0,0,0};
				if (!data.at_interp(depth_i,lat_i,lng_i,lat,lng,depth,interp)) {
					cell_out<<lat<<","<<lng<<","<<depth<<","<<interp.v<<","<<interp.u<<","<<interp.density<<"\n";
				}
			}
		}
	}
	
	return s;
}

struct ROVParameters {
	double lat,lng,depth;
	double max_horizontal, max_vertical;
	double search_radius;
	double start_time;
};

constexpr double time_penalty = 100.0/3600;

struct PointCloud {
	vector<Point> pts;
	size_t kdtree_get_point_count() const { return pts.size(); }
	double kdtree_get_pt(const size_t idx, const size_t dim) const {
		if (dim == 0) return pts[idx].lat;
		else if (dim == 1) return pts[idx].lng;
		else return pts[idx].depth;
	}
	template<class BBox>
	bool kdtree_get_bbox(BBox& bb) const {return false;}
};

using KDTree = nanoflann::KDTreeSingleIndexDynamicAdaptor<
	nanoflann::L2_Simple_Adaptor<double, PointCloud, double, int>,
	PointCloud, 3, int>;

void rov(Simulation& sim, Data& data, ROVParameters params) {
	double cur_lat = params.lat, cur_lng = params.lng, cur_depth = params.depth;
	double cur_t = params.start_time;

	list<pair<int,vector<Point>>> paths;
	for (int i=0; i<sim.points.size(); i++)
		paths.emplace_back(i, sim.points[i]);
	vector<double> time_found(sim.points.size(), inf);

	ofstream rov_path("./rov.csv");
	rov_path<<"t,lat,lng,depth\n";

	for (int iter=0; iter<100; iter++) {
		cout<<"iteration "<<iter<<", time "<<cur_t<<"\n";
		cout<<paths.size()<<" paths remaining\n";

		rov_path<<cur_t<<","<<cur_lat<<","<<cur_lng<<","<<cur_depth<<"\n";

		double best_score=-inf, arrival, bearing;
		Point where;

		auto cell = data.get_cell(cur_lat, cur_lng, cur_depth);
		if (!cell.exists) {
			cerr<<"ROV crashed\n";
			return;
		}

		int min_t_idx = ceil(cur_t/log_every_t);
		for (int t=min_t_idx; ; t++) {
			bool stop=true;
			PointCloud points;
			KDTree tree(3, points);

			vector<vector<int>> children;
			vector<int> pt_to_path;
			int path_i=-1;
			for (auto& [idx,path]: paths) {
				path_i++;
				if (path.size()<=t) continue;

				stop=false;
				pt_to_path.push_back(path_i);
				points.pts.push_back(path[t]);
				children.push_back({path_i});
			}

			if (stop) break;

			vector<nanoflann::ResultItem<int, double>> indices_dists;
			nanoflann::RadiusResultSet<double, int> results(2*params.search_radius, indices_dists);

			for (int i=0; i<path_i; i++) {
				results.clear();

				auto pt = points.pts[i];
				auto arr = {pt.lat,pt.lng,pt.depth};
				tree.findNeighbors(results, arr.begin());

				for (auto& x: indices_dists) {
					int path = pt_to_path[x.first];
					int k = children[path].size();

					if (path!=i && x.second <= (k+1.0)/k * params.search_radius) {
						children[path].push_back(i);
						tree.removePoint(x.first);
						points.pts.push_back((points.pts[x.first]*k + pt)/(k+1));
						pt_to_path.push_back(path);
						tree.addPoints(points.pts.size()-1, points.pts.size()-1);
						break;
					}
				}
			}

			for (int i=0; i<points.pts.size(); i++) {
				double at = t*log_every_t;
				auto& p = points.pts[i];

				double dist, azi1, azi2;
				data.geo.Inverse(cur_lat,cur_lng,p.lat,p.lng,dist,azi1,azi2);
				double east = sin(azi1)*dist, north = cos(azi1)*dist;

				double horiz = hypot(north/(at-cur_t) - cell.v, east/(at-cur_t) - cell.u);
				if (horiz > params.max_horizontal
					|| abs(p.depth-cur_depth)/(at-cur_t) > params.max_vertical) continue;

				// double a0=(p.north-cell->v)*(p.north-cell->v) + (p.east-cell->u)*(p.east-cell->u) - params.max_horizontal;
				// double a1=(north-v*p.north)*(p.north-cell->v) + (east - v*p.east)*(p.east-cell->u);
				// double a2=max((north-v*p.north)*(north-v*p.north) + (east - v*p.east)*(east - v*p.east), 0.001);

				// double res = a1*a1 - 4*a0*a2;
				// if (res<1e-7) continue;
				// double r1 = (-a1 - sqrt(res))/(2*a2);
				// double r2 = (-a1 + sqrt(res))/(2*a2);

				// double d1 = p.depth - cur.depth - v*p.down;
				// double d2 = p.down - params.max_vertical;
				// double d_int = abs(d1)<1e-7 ? inf : -d2/d1;
				// double d_start = d2<0 ? 0 : d_int, d_end = d2>0 ? inf : d_int;

				// double mint = max({ d_start, r1, cur_t });
				// if (mint+cur_t>(t+1)*log_every_t || mint>r2 || mint>d_end) continue;

				double sc = 1.0*children[i].size() - time_penalty*at;
				if (sc>best_score) {
					best_score=sc, where=p;
					arrival=at, bearing=azi1;
				}
			}
		}

		if (isinf(best_score)) {
			cout<<"nothing left\n";
			break;
		}

		for (auto it = paths.begin(); it!=paths.end(); it++) {
			for (int t=min_t_idx; t<it->second.size(); t++) {
				double xt = t*log_every_t;
				if (xt>arrival) break;

				auto& p = it->second[t];

				double coeff = (xt-cur_t)/(arrival-cur_t);

				double dist;
				data.geo.Inverse(cur_lat + coeff*(where.lat-cur_lat),
					cur_lng + coeff*(where.lng-cur_lng),
					p.lat,p.lng,dist);

				if (hypot(dist, cur_depth+coeff*(where.depth-cur_depth) - p.depth) <= params.search_radius) {
					time_found[it->first]=xt;
					it=paths.erase(it);
					break;
				}
			}
		}

		cur_t=arrival;
		cur_lat=where.lat, cur_lng=where.lng, cur_depth=where.depth;
	}

	ofstream out("./times.csv");
	out<<"t\n";
	for (double x: time_found) out<<x<<"\n";
}

int main(int argc, char** argv) {
	stringstream ss;
	for (int i=0; i<argc; i++) ss<<argv[i];

	json json_params = json::parse(ifstream("./parameters.json"));
	Parameters params = {
		.lat=json_params["lat"],
		.lng=json_params["lng"],
		.depth=json_params["depth"],
		.east=json_params["velocityEast"],
		.north=json_params["velocityNorth"],
		.down=json_params["velocityDown"],
		.mass=json_params["mass"],
		.volume=json_params["volume"],
		.drag=json_params["drag"],
		.area=json_params["area"]
	};

	Data data(json_params["data"]);
	auto sim = simulate(data, params);

	ROVParameters rov_params {
		.lat=json_params["rov_lat"],
		.lng=json_params["rov_lng"],
		.depth=json_params["rov_depth"],
		.max_horizontal=json_params["rov_max_horizontal"],
		.max_vertical=json_params["rov_max_vertical"],
		.search_radius=json_params["rov_search_radius"],
		.start_time=json_params["rov_start_time"]
	};

	rov(sim, data, rov_params);
}