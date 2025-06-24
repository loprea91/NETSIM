######################################################################
# This file demonstrates how to use the python interface to netsim, doing a simple scan of a parameter space 
# and plotting the results,
# using the python interface to netsim.
# Ideally we want this example script to use as few external dependencies as possible
# with default settings, it should run in ~20 minutes on a laptop
######################################################################

import os
import netsim
import numpy as np

def simulation(ge,gi):
	output_folder = f"ge= {ge} _ gi= {gi}"
	if output_folder not in [folder[:len(f'ge= {ge} _ gi= {gi}')] for folder in os.listdir("./scan_outputs/")]: # Check if simulation has already been ran

			# Create ouput directory
			os.mkdir(f"./scan_outputs/"+output_folder)
			net = netsim.netsim(params="./tests/test1.parameters", output_dir="./scan_outputs/"+output_folder,compile_on_run=False)
			net.set_params(ge=ge, gi=gi)

			#Save the parameter file
			f = open(f"./scan_outputs/"+output_folder+f"/params_ge= {ge} _ gi= {gi}.txt", 'a')

			for a in net.params:
				f.write(str(a)+" = "+str(net.params[a])+"\n")
			f.close()

			# Generate parameter file
			net.write_paramfile(net.params)

			# Run simulation
			net.run()

if __name__ == "__main__":
    
	os.system("make clean")
	os.system("make COMPILE_TYPE=performance")
	try:
		os.mkdir("./scan_outputs/")
	except:
		pass
        
	# Generate parameters sets array
	ge_values= np.linspace(0.5e-9,1e-9,2) 
	gi_values = np.linspace(1e-8,5e-8,2) 
	ge_gi_values=np.array((np.meshgrid(ge_values,gi_values))).T.reshape(-1, 2)

	#if you want it shuffled
	#np.random.shuffle(ge_gi_values)

	for ge,gi in ge_gi_values:
		simulation(ge,gi)
