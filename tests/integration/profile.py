import pyNetsim.netsim as netsim
import os
import subprocess

def profile_netsim():
    #create the object
    #PROFILE TEST1
    net = netsim.netsim(params="./tests/test1.parameters", compilekwargs=["COMPILE_TYPE=profiler"], compile_on_run=False, verbose=True)
    #net.set_params(dict(N=50000,))
    net.compile()
    #run the simulation
    os.environ["CPUPROFILE"] = "/tmp/test1.out"
    os.environ["CPUPROFILE_FREQUENCY"] = "500"
    net.run()
    #get the results
    results = net.get_results()

    net = netsim.netsim(params="./tests/test13.parameters", compilekwargs=["COMPILE_TYPE=profiler", "FAST_BINOMIAL=yes"], compile_on_run=False, verbose=True)
    #net.set_params(dict(N=50000,))
    net.compile()
    #run the simulation
    os.environ["CPUPROFILE"] = "/tmp/test13.out"
    os.environ["CPUPROFILE_FREQUENCY"] = "500"
    net.run()
    #get the results
    results = net.get_results()
    #remove the temp dir

if __name__ == "__main__":
    profile_netsim()