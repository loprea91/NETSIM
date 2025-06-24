import os
import pyNetsim.netsim as netsim
import numpy as np
import matplotlib.pyplot as plt
from pyNetsim.utils import nslfp

def test_lfp():
    #create a netsim object
    net = netsim.netsim(params='./tests/test1.parameters')
    net.run()
    #get the results
    results = net.get_results(reshape=True)
    #store the results
    params_with_matlab = netsim._add_matlab_params(net.params)
    #compute the lfp
    lfp = nslfp(params_with_matlab, results['voltages'], results['ge'], results['gi'])
    #plot the lfp
    plt.imshow(lfp.T, aspect='equal')
    plt.show()
    fig, ax = plt.subplots(2, 1)
    t = np.linspace(0, 1000, len(results['voltages'][3622, :]))
    ax[0].plot(t, results['voltages'][3622, :])
    ax[1].plot(t, results['ge'][3622,:])
    ax[0].scatter(results['spiketimes'][results['spikeids']==3622]*1000, np.ones(len(results['spiketimes'][results['spikeids']==3622]))*-0.05, color='r')
    plt.show()

if __name__ == '__main__':
    test_lfp()