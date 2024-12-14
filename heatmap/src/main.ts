import Papa from "papaparse";
import outCsv from "../../out.csv?raw";
import rovCsv from "../../rov.csv?raw";
import timesCsv from "../../times.csv?raw";
import L, { HeatLatLngTuple } from "leaflet";
import GUI from 'lil-gui';
import "leaflet.heat";
import "../node_modules/leaflet/dist/leaflet.css";

const gui = new GUI({ width: 400 });
const obj = {
  allTime: 0,
  rovTime: 0,
  time: "rov" as "all"|"rov",
  disappearPaths: true,
  blur: 5,
  radius: 4,
  max: 1,
  atTime: 0
};

const nRuns = 1000;

const lines = (outCsv as string).split("\n\n").slice(0,nRuns).map(x=>x.split("\n"));
for (let i=1; i<lines.length; i++) {
	lines[i].splice(0,0,lines[0][0]);
}

const opts: Papa.ParseConfig = {
	header: true, dynamicTyping: true,
	skipEmptyLines: true
};

const runs = lines.map(l => Papa.parse(l.join("\n"), opts).data) as {
	t: number,
	lat: number,
	lng: number,
	depth: number,
	east: number,
	north: number
}[][];

const rov = Papa.parse(rovCsv as string, opts).data as {
	t:number,lat: number,lng:number,depth:number
}[];

const times = Papa.parse(timesCsv as string, opts).data.slice(0,nRuns) as { t:number }[];

const map = L.map('map').setView([rov[0].lat, rov[0].lng], 15);

L.tileLayer.wms("https://wms.gebco.net/mapserv?", {
  // attribution: '&copy; <a href="https://osm.org/copyright">OpenStreetMap</a> contributors',
  layers: "GEBCO_LATEST"
}).addTo(map);

L.control.scale().addTo(map);

let undo: null|(()=>void) = null;
function render() {
  if (undo!=null) undo();
  const toRemove: {remove:()=>void}[]=[];
  const t = obj.time=="all" ? obj.allTime : obj.rovTime;

  const heat = L.heatLayer(runs
    .filter((_,i) => {
      return !obj.disappearPaths || times[i].t >= t;
    })
    .flatMap(x=>x.map(y=>[
      y.lat,y.lng, 50.0/(1.0 + (obj.atTime/1000)*(y.t-t)*(y.t-t))
    ] satisfies HeatLatLngTuple)), obj).addTo(map);
  toRemove.push(heat);

  const rovPts: [number,number][] = [];
  for (let i=0; i<rov.length; i++) {
    if (rov[i].t>t) {
      if (i==0) break;

      const c = (t-rov[i-1].t)/(rov[i].t-rov[i-1].t);
      rovPts.push([
        c*rov[i].lat + (1-c)*rov[i-1].lat,
        c*rov[i].lng + (1-c)*rov[i-1].lng
      ]);

      break;
    } else {
      rovPts.push([rov[i].lat,rov[i].lng]);
    }
  }

  if (rovPts.length>0) {
    const line = L.polyline(rovPts, {
      color: "red",
      weight: 5
    }).addTo(map);

    const circ = L.circle(rovPts[rovPts.length-1], {
      radius: 200, fillColor: "red",
      fillOpacity: 0.3,
      color: "red", weight: 2
    }).addTo(map);

    toRemove.push(line,circ);
  }

  undo=() => {
    for (const x of toRemove) x.remove();
  };
}

gui.onChange(render);
gui.add(obj, "rovTime", rov[0].t, rov[rov.length-1].t);

gui.add(obj, "allTime", 0, runs.flatMap(r=>r.map(x=>x.t)).reduce((a,b)=>Math.max(a,b)));
gui.add(obj, "atTime", 0, 1);
gui.add(obj, "time", ["all", "rov"]);

gui.add(obj, "blur", 0, 10);
gui.add(obj, "radius", 0, 10);
gui.add(obj, "max", 0, 10);
gui.add(obj, "disappearPaths");

render();
