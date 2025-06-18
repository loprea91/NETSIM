import numpy as np
import os
import pyNetsim.netsim as netsim
import time

import matplotlib.pyplot as plt
def test1():
    #create the object
    net = netsim.netsim(params="/root/NETSIM-dev/tests/test13.parameters", compilekwargs=["COMPILE_TYPE=performance"], compile_on_run=False, verbose=True)
    net.compile()
    #sample some output codes 
    codes = np.random.randint(0, 100, 5)
    old_binom = []
    old_binom_times = []
    ge = []
    #run the simulation
    for c in codes:
        st = time.time()
        net.set_params(output_file_code=str(c).zfill(3))
        net.run()
        #get the results
        results = net.get_results(reshape=True)
        #check the spike amount
        spike_count = len(results["spikesids"])
        #make a eventplot
        plt.eventplot(results['spiketimes'], lineoffsets=spike_count, linelengths=0.5)
        plt.show()
        old_binom.append(spike_count)
        old_binom_times.append(time.time()-st)
    
    #now run the simulation with the fast binomial
    net = netsim.netsim(params="/root/NETSIM-dev/tests/test13.parameters", compilekwargs=["COMPILE_TYPE=performance", "FAST_BINOMIAL=yes"], compile_on_run=False, verbose=True)
    net.compile()
    fast_binom = []
    fast_binom_times = []
    for c in codes:
        st = time.time()
        net.set_params(output_file_code=str(c).zfill(3))
        net.run()
        #get the results
        results = net.get_results()
        #check the spike amount
        spike_count = len(results["spikesids"])
        fast_binom.append(spike_count)
        fast_binom_times.append(time.time()-st)
    #make a dataframe of the results
    import pandas as pd
    df = pd.DataFrame({"old_binom":old_binom, "fast_binom":fast_binom, "old_binom_times":old_binom_times, "fast_binom_times":fast_binom_times})
    df.to_csv("test1.csv")
    #plot the results
    plt.plot(old_binom, label="old binomial")
    plt.plot(fast_binom, label="fast binomial")
    plt.legend()
    plt.show()

def test2():
    #test the increasing N
    #create the object
    net = netsim.netsim(params="/root/NETSIM-dev/tests/test13.parameters", compilekwargs=["COMPILE_TYPE=performance"], compile_on_run=False, verbose=True)
    net.compile()
    #create a logspace of N
    Ns = np.logspace(3, 4, 5)
    old_binom = []
    old_binom_times = []
    #run the simulation
    for N in Ns:
        st = time.time()
        net.set_params(N=int(N))
        net.run()
        #get the results
        results = net.get_results()
        #check the spike amount
        spike_count = len(results["spikesids"])
        old_binom.append(spike_count)
        old_binom_times.append(time.time()-st)

    #now run the simulation with the fast binomial
    net = netsim.netsim(params="/root/NETSIM-dev/tests/test13.parameters", compilekwargs=["COMPILE_TYPE=performance", "FAST_BINOMIAL=yes"], compile_on_run=False, verbose=True)
    net.compile()
    fast_binom = []
    fast_binom_times = []
    for N in Ns:
        st = time.time()
        net.set_params(N=int(N))
        net.run()
        #get the results
        results = net.get_results()
        #check the spike amount
        spike_count = len(results["spikesids"])
        fast_binom.append(spike_count)
        fast_binom_times.append(time.time()-st)

    #make a dataframe of the results
    import pandas as pd
    df = pd.DataFrame({"old_binom":old_binom, "fast_binom":fast_binom, "old_binom_times":old_binom_times, "fast_binom_times":fast_binom_times})
    df.to_csv("test2.csv")
    #plot the results
    plt.plot(old_binom, label="old binomial")
    plt.plot(fast_binom, label="fast binomial")
    plt.legend()


if __name__ == "__main__":
    test1()