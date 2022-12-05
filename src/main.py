
def main(topo, traj, pore_params, grid_params, box_params, cutoff, volt_params, param_path, out, log_flag):


	from sem_params import SEMParams
#	from conductivity import ConductivityOf
	from conductivity import Residue2Conductivity, set_conductivity_parameters
	from molecules import Molecules
	from sem import StericExclusionModel
	from const import SEMConst
	import numpy as np
	import time
	from tqdm import tqdm

	start_time = time.time()

	SEMParams.grid_params  = np.array(grid_params)
	SEMParams.pore_params  = np.array(pore_params)
	SEMParams.box_params   = np.array(box_params)
	SEMParams.cutoff = cutoff
	SEMParams.volt_bottom  = volt_params[0]
	SEMParams.volt_top = volt_params[1]
	SEMParams.conductivity_path = param_path
	SEMParams.log_flag = log_flag
	params = SEMParams

	
	#conductivity_of = ConductivityOf
	#ConductivityOf.set_parameter(param_path)
	set_conductivity_parameters()

	time_series = Molecules(topo, traj)

	model = StericExclusionModel(params)
	#model = StericExclusionModel(params, conductivity_of)


	if SEMParams.log_flag:
		print("Computing electric current of each frame ...")

	currents = list()
	step_iterator = tqdm(range(time_series.framenum)) if SEMParams.log_flag else range(time_series.framenum)
	for istep in step_iterator:
	#for istep in range(250, 300):
		positions = time_series.frame(istep)
		resnames  = time_series.resnames
		flow_top, flow_btm = model.run(positions, resnames)
		currents.append([flow_top, flow_btm])
		#if SEMParams.log_flag and (istep % SEMConst.LOG_DUMP_STEP == 0):
		#	computation_time = time.time()
		#	print("{:8d} Steps - Elapsed time: {:16.8} sec.".format(istep, computation_time - start_time))

	if SEMParams.log_flag:
		print("Done")
	
	with open(out, 'w') as ofs:
		ofs.write("#     Step      Current_Top   Current_Bottom\n")
		for istep, val in enumerate(currents):
			ofs.write("{:10d} {:16.8e} {:16.8e}\n".format(istep, val[0], val[1]))


if __name__ == "__main__":

	import argparse
	import toml

	parser = argparse.ArgumentParser()
	parser.add_argument("--toml", required = True)
	parser.add_argument("--out", default = "current.dat")
	parser.add_argument("-v", action = "store_true")
	args = parser.parse_args()
	inp = toml.load(open(args.toml, 'r'))

	fnames = inp["filename"]
	pore_params = [inp["parameter"]["pore"]["length"], inp["parameter"]["pore"]["upper"], inp["parameter"]["pore"]["lower"]]
	grid_params = [inp["parameter"]["grid"]["x"], inp["parameter"]["grid"]["y"], inp["parameter"]["grid"]["z"]]
	box_params  = [inp["parameter"]["box"]["x"], inp["parameter"]["box"]["y"], inp["parameter"]["box"]["z"]]
	cutoff = inp["parameter"]["cutoff"]
	volt_params = [inp["parameter"]["voltage"]["bottom"], inp["parameter"]["voltage"]["top"]]

	main(fnames["topo"], fnames["traj"], pore_params, grid_params, box_params, cutoff, volt_params, fnames["conductivity"], args.out, args.v)
