########
# This file contains functions for loading and processing data from the simulation
# primairly we are porting the matlab functions from the analysis folder
# We are also adding some functions for plotting and analysis
########

from scipy.signal import welch, vectorstrength
import scipy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
import copy
import math

def zscore(x, axis):
    return (x - np.mean(x, axis=axis, keepdims=True)) / np.std(x, axis=axis, keepdims=True)

def nslfp(params, v, ge, gi):
    v = v[0:int(params['nbins_exc']), :]
    ge = ge[0:int(params['nbins_exc']), :]
    gi = gi[0:int(params['nbins_exc']), :]
    tau = 0.006 / (params['dt'] * params['record_downsample_factor'])
    exc = np.roll(ge * (params['Ee'] - v), int(tau), 1)
    inh = gi * (params['Ei'] - v)
    lz = zscore(abs(exc - 1.65 * inh), 1)
    return lz

def binarize_spike_train(spike_train, dt=1e-4, run_length=1):
    """_summary_

    Args:
        spike_train (_type_): _description_
        dt (_type_): _description_
        run_length (_type_): _description_

    Returns:
        _type_: _description_
    """
    #get the number of bins
    n_bins = int(run_length/dt)
    #create the bins
    bins = np.arange(0, run_length, dt)
    #bin the spike train
    binned_spike_train = np.histogram(spike_train, bins=bins)[0]
    return binned_spike_train