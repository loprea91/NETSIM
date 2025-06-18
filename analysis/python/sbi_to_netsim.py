import os
import pyNetsim.netsim as netsim
from pyNetsim.utils import binarize_spike_train
import numpy as np
import torch
import matplotlib.pyplot as plt
import pickle as pkl
import pyNetsim.utils as netsimutils

import scipy.signal as signal
import scipy.fft as fft
#sbi imports
from sbi import utils as utils
from sbi import analysis as analysis
#import SNPE from sbi
from sbi.inference import SNPE
from sbi.inference.base import infer, simulate_for_sbi, prepare_for_sbi
import subprocess
from functools import partial

def get_isi(spike_times):
    #compute the diff along the last axis
    diff = [np.diff(x, axis=-1) for x in spike_times if isinstance(x, np.ndarray)]
    return diff


def run_and_plot(params):
    #run the simulation with the most likely params
    net = netsim.netsim()
    param_dict = {'ge': params[0], 'gi': params[1], 'max_conduction_speed': params[2]} #create a dictionary of the params
    net.set_params(param_dict) #set the params
    net.run() #run the simulation
    results = net.get_results(reshape=True, reshape_spikes=True) #get the results, returns a dictionary of lists/ndarrays
    #plot the results
    #spike raster
    plt.figure()
    plt.eventplot(results['spiketimes'])
    
    log_isi = np.hstack(get_isi(results['spiketimes']))*1000
    plt.figure()
    plt.hist(log_isi, bins=np.logspace(0, 3, 100))
    


def error_function(net, params):
    outputdir = './temp' + str(os.getpid()) + '/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    net.outputdir = outputdir
    #turn the params into a dictionary
    params = params.numpy()
    #empty lists to hold the results
    obs_out =[]
    #if the params are 2 dimensional, then we are running a single simulation
    batch_size = params.shape[0] if len(params.shape) > 1 else 1
    for i in range(batch_size):
        if len(params.shape)>1: #if we are running a batch of simulations
            param_dict = {'ge': params[i][0], 'gi': params[i][1], 'max_conduction_speed':params[i][2]} #create a dictionary of the params
        else:
            param_dict = {'ge': params[0], 'gi': params[1], 'max_conduction_speed':params[2]} #create a dictionary of the params
        net.set_params(param_dict) #set the params
        try:
            net.run() #run the simulation
            results = net.get_results(reshape=True, reshape_spikes=True) #get the results, returns a dictionary of lists/ndarrays
        except:
            #make some zeros if the simulation failed
            results = {'spiketimes': np.zeros((net.params['N'])), 'voltage': np.zeros((net.params['N']))}
        
        n_neurons = net.params['N'] #get the number of neurons
        isi = get_isi(results['spiketimes']) #get the isi
        if len(isi) == 0:
            #if there are no spikes, then make some zeros
            hist = np.zeros(99) + 1e-15
        else:
            log_isi = np.hstack(isi)*1000
            #compute the histogram of the log isi
            hist, bin_edges = np.histogram(log_isi, bins=np.logspace(0, 3, 100))
        obs_out.append(hist)
    
    
    #normalize the log isi histogram
    obs_out = np.array(obs_out)
    obs_out = obs_out / np.sum(obs_out, axis=-1, keepdims=True)
    #delete the temp dir
    subprocess.run(["rm", "-r", outputdir])
    return torch.tensor(obs_out, dtype=torch.float32) #return the results

def main(load_previous=True):
    #init the netsim object
    net = netsim.netsim(compile_on_run=False, verbose=True)
    net.set_params({'N': 1000, 'K':100, 'ge': 0.1, 'gi': 0.1, 'max_conduction_speed': 0.1, 'poisson_rate': 10}) #set the number of neurons
    net.compile() #compile the simulation
    #create the function that will be used by sbi
    error_function_part = partial(error_function, net)
    if load_previous:
        #load the previous density estimator with pickle
        with open('density_estimator.pkl', 'rb') as f:
            density_estimator = pkl.load(f)
        with open('posterior.pkl', 'rb') as f:
            posterior = pkl.load(f)
        #create a sample by drawing from a normal distribution with mean 100 and std 10
        sample = np.hstack((np.random.normal(100, 10, 1000), np.random.normal(100, 10, 1000)))
        #histogram of the log isi
        hist, _ = np.histogram(sample, bins=np.logspace(0, 3, 100))
        #normalize the histogram
        hist = hist / np.sum(hist)
        #convert to a tensor
        observation = torch.tensor(hist, dtype=torch.float32)
        samples = posterior.sample((10000,), x=observation) #sample from the posterior
        log_probability = posterior.log_prob(samples, x=observation) #get the log probability of the samples
        _ = analysis.pairplot(samples, figsize=(6, 6), labels=['ge', 'gi', 'delay'], limits=[[0.1e-11,4e-9], [1e-10, 30e-9], [0.02, .1]],) #plot the posterior
        plt.show()
    else:
        prior = utils.BoxUniform(low=[0.1e-11, 1e-10, 0.02], high=[4e-9, 30e-9, .1]) #define the prior, first param is ge, second param is gi

        simulator, prior = prepare_for_sbi(error_function_part, prior) #prepare the simulator for sbi, some preprocessing
        #simulate data
        num_simulations = 10000
        theta, x = simulate_for_sbi(simulator, prior, num_simulations, simulation_batch_size=2, num_workers=8)
        #run inference
        inference = SNPE(prior)
        density_estimator = inference.append_simulations(theta, x).train() #train the density estimator
        posterior = inference.build_posterior(density_estimator) #build the posterior
        #save the density estimator and posterior
        with open('density_estimator.pkl', 'wb') as f:
            pkl.dump(density_estimator, f)
        with open('posterior.pkl', 'wb') as f:
            pkl.dump(posterior, f)
        #plot the posterior
        observation = x[5] #define the observation
        samples = posterior.sample((10000,), x=observation) #sample from the posterior
        log_probability = posterior.log_prob(samples, x=observation) #get the log probability of the samples
        _ = analysis.pairplot(samples, figsize=(6, 6)) #plot the posterior
    plt.savefig("/root/NETSIM-dev/posterior.png")
    #grab the most likely parameters
    most_likely_params = posterior.set_default_x(observation).map()
    
    run_and_plot(most_likely_params) #run the simulation with the most likely params and plot the results

    #also plot the least likely params
    least_likely_params = samples[np.argmin(log_probability)]
    run_and_plot(least_likely_params) #run the simulation with the least likely params and plot the results

    #also plot the median params
    median_params = samples[np.argmin(np.abs(log_probability - np.median(log_probability)))]
    run_and_plot(median_params) #run the simulation with the median params and plot the results

    #also plot for a 400 isi value
    #sample = np.hstack((np.random.normal(20, 10, 1000), np.random.normal(200, 10, 1000)))
    #hist, _ = np.histogram(sample, bins=np.logspace(0, 3, 100))
    #hist = hist / np.sum(hist)
    #observation = torch.tensor(hist, dtype=torch.float32)
    #samples = posterior.sample((10000,), x=observation) #sample from the posterior
    #log_probability = posterior.log_prob(samples, x=observation) #get the log probability of the samples
    #most_likely_params = posterior.set_default_x(observation).map()
    #run_and_plot(most_likely_params) #run the simulation with the median params and plot the results
    plt.hist(sample, bins=np.logspace(0, 3, 100), color='r', alpha=0.5, density=True)
    plt.show()
    

    

if __name__ == '__main__':
    main()
