import os
import pyNetsim.netsim as netsim
import numpy as np

def main():
    #create a netsim object
    net = netsim.netsim(verbose=True)

    #run the simulation
    net.run()
    #get the results
    net.get_results(reshape=True)
    #export results to matlab mat files
    net.export_to_mat("./results")

if __name__ == '__main__':
    main()
