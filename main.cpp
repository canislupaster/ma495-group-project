#include <cmath>
#include <highfive/H5File.hpp>
#include <initializer_list>
#include <limits>
#include <random>
#include <set>
#include <map>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <nlohmann/json.hpp>

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

	HighFive::File file(json_params["data"]);
	auto dataset = file.getDataSet("data");
	auto data = dataset.read<vector<array<double,unsigned(Cols::num_column)>>>();
	map<double, int, EpsCmp> lat_idx, lng_idx, depth_idx;

	for (auto& v: data) {
		depth_idx.insert({v[int(Cols::depth)],-1});
		lat_idx.insert({v[int(Cols::latitude)],-1});
		lng_idx.insert({v[int(Cols::longitude)],-1});
	}

	vector<double> depths, lats, lngs;
	for (auto [a,_]: depth_idx) depths.push_back(a);
	for (auto [a,_]: lat_idx) lats.push_back(a);
	for (auto [a,_]: lng_idx) lngs.push_back(a);

	for (auto map: {&lat_idx, &lng_idx, &depth_idx}) {
		int i=0;
		for (auto& [a,b]: *map) b=i++;

		map->insert({-inf, -1});
		map->insert({inf, -1});
	}

	int num_lat=lats.size(), num_lng=lngs.size(), num_depth=depths.size();

	vector<Cell> cells(num_lat*num_lng*num_depth);
	auto at = [&cells,num_lat,num_lng](int lat_i, int lng_i, int depth_i) -> Cell& {
		return cells[lat_i + num_lat*(lng_i + num_lng*depth_i)];
	};

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

	constexpr double dt = 5;
	constexpr double maxt = 24*3600*5; //a day
	constexpr int log_every = 100;
	constexpr double noise_update_second = 0.001;
	constexpr double g = 9.81;

	ofstream out("./out.csv");

	auto const& geo = Geodesic::WGS84();

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

	auto at_interp = [&](int depth_i, int lat_i, int lng_i, double cur_lat, double cur_lng, double cur_depth, Cell& interpolated) -> bool {
		double tot = (depths[depth_i]-depths[depth_i-1])
			*(lats[lat_i]-lats[lat_i-1])
			*(lngs[lng_i]-lngs[lng_i-1]);

		for (int z=0; z<=1; z++) {
			for (int y=0; y<=1; y++) {
				for (int x=0; x<=1; x++) {
					auto& cell = at(lat_i+x, lng_i+y, depth_i+z);
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


	BBox box;

	for (int sim_i=0; sim_i<1000; sim_i++) {
		double t=0;
		int i=0;
		bool crashed=false, exited=false;

		double cur_lat = params.lat, cur_lng = params.lng, cur_depth = params.depth;
		double cur_down = params.down, cur_east = params.east, cur_north = params.north;

		double down_noise=down_dist(rng), u_noise=u_dist(rng), v_noise=v_dist(rng);

		int depth_i = depth_idx.upper_bound(cur_depth)->second;
		int lat_i = lat_idx.upper_bound(cur_lat)->second;
		int lng_i = lng_idx.upper_bound(cur_lng)->second;

		while (t<maxt) {
			if ((i++)%log_every == 0) {
				initializer_list<double> props = {t,cur_lat,cur_lng,cur_depth,cur_east,cur_north};
				for (int prop_i=0; prop_i<props.size(); prop_i++) {
					out<<*(props.begin()+prop_i)<<",\n"[prop_i==props.size()-1];
				}

				box.add({cur_lat,cur_lng,cur_depth});
			}

			t+=dt;

			cur_depth += dt*cur_down;

			double ang = atan2(cur_east, cur_north)*180/M_PI;
			if (!isnan(ang)) {
				if (ang<0) ang=360-ang;
				geo.Direct(cur_lat, cur_lng, ang, dt*hypot(cur_east, cur_north), cur_lat, cur_lng);
			}

			if (iter(depth_i, depths, cur_depth)) {
				crashed=true; break;
			}
			if (iter(lat_i, lats, cur_lat)
				|| iter(lng_i, lngs, cur_lng)) {
				exited=true; break;
			}

			double a = noise_update_second*dt, b=1-a;
			down_noise = a*down_dist(rng) + b*down_noise;
			u_noise = a*u_dist(rng) + b*u_noise;
			v_noise = a*v_dist(rng) + b*v_noise;

			Cell interpolated = {.density=0, .u=u_noise, .v=v_noise};
			if (at_interp(depth_i,lat_i,lng_i,cur_lat,cur_lng,cur_depth,interpolated)) {
				crashed=true;
				break;
			}

			double mul = dt/params.mass;

			cur_down += mul*(
				g*(params.mass - interpolated.density*params.volume)
				+ sgn(down_noise - cur_down)*interpolated.density*(down_noise - cur_down)*(down_noise - cur_down)*params.drag*params.area/2
			);

			cur_north += mul*(
				sgn(interpolated.v - cur_north)*interpolated.density*(interpolated.v - cur_north)*(interpolated.v-cur_north)*params.drag*params.area/2
			);

			cur_east += mul*(
				sgn(interpolated.u - cur_east)*interpolated.density*(interpolated.u - cur_east)*(interpolated.u - cur_east)*params.drag*params.area/2
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

	for (double depth=max(box.b_min[2], depths[0]); depth<=min(box.b_max[2], depths.back()); depth+=depth_inc) {
		iter(depth_i, depths, depth);
		for (double lat=max(box.b_min[0], lats[0]); lat<=min(box.b_max[0], lats.back()); lat+=lat_inc) {
			iter(lat_i, lats, lat);
			for (double lng=max(box.b_min[1], lngs[0]); lng<=min(box.b_max[1], lngs.back()); lng+=lng_inc) {
				iter(lng_i, lngs, lng);
				Cell interp {0,0,0,0};
				if (!at_interp(depth_i,lat_i,lng_i,lat,lng,depth,interp)) {
					cell_out<<lat<<","<<lng<<","<<depth<<","<<interp.v<<","<<interp.u<<","<<interp.density<<"\n";
				}
			}
		}
	}
}