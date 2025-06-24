#!/usr/bin/python3

import sys
from pathlib import Path
import subprocess
import os
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

try:
    import joblib
    PARALELL = True
except:
    print("joblib not found, test for parallelization will be skipped")
    PARALELL = False

if len(sys.argv) > 1:
    print("Using " + sys.argv[1] + " is the root NETSIM directory\n")
    ROOT_PATH = Path(sys.argv[1]).absolute()
else:
    print("Assuming the cwd is the root NETSIM directory\n")
    ROOT_PATH = Path.cwd().absolute()
    
#sys.path.append(str(Path(ROOT_PATH) / "python_interface/"))
import pyNetsim.netsim as netsim

def test_pathing():
    ################## test calling a netsim object from the same dir, outputting to a different dir ##################
    #here we are testing the output dir is created and the results are saved to it
    # the output dir is a temp dir based on the process id
    # ideally we want that the python interface is telling the C code where to save the results
    # and not the C code saving the results to the specified output dir and not the cwd
    #path to NETSIM dir
    print("Testing netsim pathing ", end="")
    root_path = str(ROOT_PATH)+"/"

    #create a temp dir to run the test in
    outputdir = './temp' + str(os.getpid())
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    #make a netsim obj
    net = netsim.netsim(params=root_path+"tests/test1.parameters",
                        output_dir=outputdir,
                        netsim_dir=root_path,
                        compile_on_run=True,
                        verbose=False)
    #run the simulation
    net.run()
    #enforce the resulst are saved to the output dir
    assert(os.path.exists(outputdir+"/00000002ge.bin"))
        
    #get the results
    results = net.get_results()

    #remove the temp dir
    subprocess.run(["rm", "-rf", outputdir])

    ################## test calling a netsim object from a different dir ##################
    root_path = str(ROOT_PATH)+"/"
    #create a temp dir to run the test in
    outputdir = './temp' + str(os.getpid()) + '/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    #copy src to the temp dir
    subprocess.run(["cp", "-r", root_path+"src", outputdir])
    #copy python_interface to the temp dir
    subprocess.run(["cp", "-r", root_path+"python_interface", outputdir])
    #copy tests to the temp dir
    subprocess.run(["cp", "-r", root_path+"tests", outputdir])
    #copy makefile to the temp dir
    subprocess.run(["cp", root_path+"Makefile", outputdir])

    #make a netsim obj
    net = netsim.netsim(params=outputdir+"/tests/test1.parameters",
                        output_dir=outputdir,
                        netsim_dir=outputdir,
                        compile_on_run=True,
                        verbose=False)

    #run the simulation
    net.run()
    #enforce the resulst are saved to the output dir
    assert(os.path.exists(outputdir+"/00000002ge.bin"))

    #get the results
    results = net.get_results(reshape=True)

    #remove the temp dir
    subprocess.run(["rm", "-rf", outputdir])

    ####### finally test calling from a different dir ######

    working_dir = os.getcwd()
    os.chdir(os.path.expanduser("~"))
    #make a netsim obj
    net = netsim.netsim(params=root_path+"/tests/test1.parameters",
                        netsim_dir=root_path,
                        compile_on_run=True,
                        verbose=False)

    net.set_params(output_file_code="003")

    #run the simulation
    net.run()
    #enforce the results are saved to the output dir
    assert(os.path.exists("00000003ge.bin"))
    print("\u001b[32m" + "Successful" + "\u001b[0m")
    
    #change back to the working dir
    os.chdir(working_dir)
    return True
    
def test_parallel():
    print("Testing netsim parallelism ", end="")
        
    def internal_call(x=[1e-11, 1e-9]):
        root_path = str(ROOT_PATH)+"/"
        #create a folder for the PID
        outputdir = './temp' + str(os.getpid()) + '/'
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        #copy src to the temp dir
        subprocess.run(["cp", "-r", root_path+"src", outputdir])
        #copy python_interface to the temp dir
        subprocess.run(["cp", "-r", root_path+"python_interface", outputdir])
        #copy tests to the temp dir
        subprocess.run(["cp", "-r", root_path+"tests", outputdir])
        #copy makefile to the temp dir
        subprocess.run(["cp", root_path+"Makefile", outputdir])

        #make a netsim obj
        net = netsim.netsim(params=outputdir+"/tests/test1.parameters",
                            output_dir=outputdir,
                            netsim_dir=outputdir,
                            compile_on_run=True,
                            verbose=False)

        net.set_params(ge=x[0], gi=x[1])
        #run the simulation
        net.run()
        
        get_results = net.get_results(reshape=False)
        #remove the temp dir
        subprocess.run(["rm", "-rf", outputdir])
        return len(get_results['spikeids'])

    #test parallel
    import numpy as np
    from joblib import Parallel, delayed
    #create a sample parameter space of ge values
    ge_values = np.linspace(1e-11, 4e-11, 5) # in siemens
    gi_values = np.linspace(1e-9, 8e-9, 5) # in ms
    #create a array of all possible combinations of ge and gi
    ge_grid, gi_grid = np.meshgrid(ge_values, gi_values)
    ge_gi_values = np.array((ge_grid, gi_grid)).T.reshape(-1, 2) #cheap hack to get the right shape,
    #run the simulation in parallel, using joblib
    results = Parallel(n_jobs=3)(delayed(internal_call)(x) for x in ge_gi_values)
    #check the results
    #make sure all the results are not the same
    if(len(results) == len(ge_gi_values)) or (len(set(results)) > 1):
        print("\u001b[32m" + "Successful" + "\u001b[0m")
    else:
        print("\u001b[31m" + "Failed" + "\u001b[0m")
        
def test_results():
    #path to NETSIM dir
    root_path = str(ROOT_PATH)+"/"
    
    inputs = ("test1.parameters",    #Random connect
              "test2.parameters",    #Gaussian connect
              #"test3.parameters",    #Test connect
              "test4.parameters",    #Gaussian connect
              "test5.parameters",    #Gaussian connect with gamma distribution
              "test6.parameters",    #Gaussian connect
              "test7.parameters",    #Gaussian connect 2D
              #"test8.parameters",    #Parameter scans
              #"test9.parameters",    #Parameter scans
              #"test10.parameters",   #Multilayer connect
              #"test10-2.parameters", #Multilayer connect
              "test11.parameters",   #Current Injection
              #"test12.parameters",   #Clustered connect
              "test13.parameters",   #Poisson inputs
              "test14.parameters",   #Variable Poisson inputs
              "test15.parameters",   #Sparse current injection
              "test16.parameters",   #long simulation
              )

    job = 0
    i = 0

    print("Running Test cases")
    for test in inputs:
        print("Testing netsim with parameter file: " + test + " - Result: ", end="")
        net = netsim.netsim(params=root_path+"tests/" + test,
                            output_dir="./data",
                            netsim_dir=root_path,
                            compile_on_run=True,
                            verbose=False)
        net.run()
        results = net.get_results()
        produced_spikes = len(results["spikeids"])
        
        #get true spike count    
        f = open(root_path + "tests/" + test)
        for line in f:
            if line.find("total spikes:") != -1:
                expected_spikes = line.strip(" #\n").split(" ")[2]
        f.close()

        try:
            expected_spikes = int(expected_spikes)
            produced_spikes = int(produced_spikes)
            if produced_spikes == expected_spikes:
                print("\u001b[32m" + "Successful" + "\u001b[0m")
            else:
                print("\u001b[31m" + "Failed" + "\u001b[0m")
        except:
            print("\u001b[31m" + "Failed" + "\u001b[0m")
            
def main():
    #test_pathing()
        
    #if PARALELL:
    #    test_parallel()

    test_results()

if __name__ == "__main__":
    main()
