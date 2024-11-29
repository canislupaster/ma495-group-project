import os
import copernicusmarine
from h5py import File
from dotenv import load_dotenv
import functools

load_dotenv()

dataset_to_vars = {
	"cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m": ["so"], # practical salinity
	"cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m": ["thetao"], # potential temperature,
	"cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m": ["uo","vo"] # uv
}

with File("./data.hdf5", "w") as hdf:
	print("loading datasets")

	first=True
	frames=[]

	for dataset, vs in dataset_to_vars.items():
		print(f"loading variables {", ".join(vs)} from {dataset}")

		frame = copernicusmarine.read_dataframe(
			dataset_id=dataset,
			variables=vs,
			minimum_longitude=11.271825087079575,
			maximum_longitude=25.50476714555446,
			minimum_latitude=34.95555342377165,
			maximum_latitude=41.58203792312717,
			start_datetime="2024-12-07T00:00:00",
			end_datetime="2024-12-07T00:00:00",
			username=os.environ["user"],
			password=os.environ["pass"]
		)

		print(f"loaded {dataset}")

		frame=frame.reset_index()[(["latitude", "longitude", "depth"] if first else []) + vs]
		first=False

		frames.append(frame)

	# assume all the indices are the same, which seems to be the case :)
	hdf.create_dataset("data", data=functools.reduce(lambda x,y: x.join(y), frames)[["depth", "latitude", "longitude", "so", "thetao", "uo", "vo"]])
