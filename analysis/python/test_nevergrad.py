import os
import pyNetsim.netsim as netsim
import numpy as np
import nevergrad as ng
import scipy.signal as signal
from joblib import Parallel, delayed
from pyNetsim.utils import nslfp
from matplotlib import pyplot as plt

def error_func(netsim_obj, params):
    """This is the error function that we will use to optimize the parameters
        takes in a netsim object and a dictionary of parameters
        returns a single value, the error"""
    #since we are using the ask and tell interface, we need to pass the netsim object
    #and the params separately
    #set the params
    param_dict = {'ge': params['ge'].value * 1e-9 # to convert from nS to S
                  , 'gi': params['gi'].value * 1e-9 # to convert from nS to S
                  , 'max_conduction_speed': params['max_conduction_speed'].value, 
                  'poisson_rate': params['poisson_rate'].value,
                  'connector': params['connector'].value}
    netsim_obj.set_params(param_dict)
    #create a dir for the output based on param space
    outputdir = os.path.expanduser('~') + '/temp/' + str(os.getpid()) + '/'
    netsim_obj.outputdir = outputdir
    #run the simulation
    r = netsim_obj.run()
    if r < 0:
        return 999
    #get the results
    results = netsim_obj.get_results(reshape=True, reshape_spikes=True)

    #### Error calculation ####
    # below you would customize to your needs, here we will just use a simple error function
    #add the matlab params
    params_with_matlab = netsim._add_matlab_params(netsim_obj.params)
    #compute the lfp
    lfp = nslfp(params_with_matlab, results['voltages'], results['ge'], results['gi']) 
    #remove the first 200 ms
    lfp = lfp[0][int(0.200/(params_with_matlab['dt']*params_with_matlab['record_downsample_factor'])):]
    #compute the power spectrum
    f, Pxx = signal.welch(lfp, fs=1/params_with_matlab['dt'], nperseg=512, nfft=2048)
    #find the peak frequency
    peak_freq = f[np.argmax(Pxx)]
    #we want a peak frequency of 30 Hz
    freq_error = np.abs(peak_freq - 60)
    #and an overall spike rate of 15 Hz
    spike_rate = np.sum([len(x) for x in results['spiketimes']])/params_with_matlab['N']
    rate_error = np.abs(spike_rate - 15)
    #we want the power of the lfp to to be 1
    power_error = np.abs(np.amax(Pxx) - 1)
    #in our case we will delete the output dir, but you could also keep it
    os.system('rm -rf ' + outputdir)
    print(f"rate error {rate_error}, power error {power_error}, freq error {freq_error}")
    return rate_error + freq_error + 10**power_error #return a single value, the error
    #here we are just adding the errors, but you could also use a weighted sum, or a different function
    #rate error and freq error are in Hz, power error is in V^2/Hz, so we need to scale it up
    

def main():
    #init the netsim object
    netsim_obj = netsim.netsim(compile_on_run=False, verbose=True)
    netsim_obj.set_params({'N': 1000, 'K':100, 'ge': 0.1, 'gi': 0.1, 'max_conduction_speed': 0.1, 'poisson_rate': 10})
    #compile it, this will create the executable, which will be used by the optimizer
    netsim_obj.compile()
    #build a search space in nevergrad
    #ge, gi, max_conduction_speed, poisson_rate
    params = ng.p.Dict(ge=ng.p.Scalar(lower=0.1, upper=10), #nevergrad allows you to set bounds on the parameters, also init values
                        gi=ng.p.Scalar(lower=0.5, upper=80), 
                        max_conduction_speed=ng.p.Scalar(lower=0, upper=1), 
                        poisson_rate=ng.p.Log(lower=0.5, upper=20), #here we can also set the scale at which the parameters are sampled, log scale will optimize this parameter in log space
                        connector=ng.p.Choice([1, 2]) #we can also optimize categorical parameters, here we will optimize the connector type which is either 1 or 2
                                                    #this is useful if you want to optimize a parameter that is not continuous or a string
                                                    #ng.p.TransitionChoice is also useful if you want to optimize a parameter that is not continuous or a string
                                                    # but you want to allow for a transition between values in a specific order
                        ) 
    
    #here we will use the portfolio optimizer
    optimizer = ng.optimizers.Portfolio(params, #the search space
                                        budget=400,#the number of iterations
                                        num_workers=8, #the number of workers to use
                                        ) 
    #run the optimization, here we will use the ask and tell interface, which is the most flexible
    #we could also use the minimize function, which is a one liner
    # x = optimizer.minimize(error_func, verbosity=2)
    for i in range(optimizer.budget//optimizer.num_workers):
        x_list = []
        for j in range(optimizer.num_workers): #this is the number of workers
            x = optimizer.ask()
            x_list.append(x)
        #joblibs syntax is a little weird, but this will run the error function in parallel
        y = Parallel(n_jobs=8)(delayed(error_func)(netsim_obj, x) for x in x_list)
        for x, y in zip(x_list, y):
            optimizer.tell(x, y)
        print(f"iteration {i}, best value {optimizer.provide_recommendation()}")
    #get the best params
    best_params = optimizer.provide_recommendation()
    #set the params
    param_dict = {'ge': best_params['ge'].value * 1e-9 # to convert from nS to S
                    , 'gi': best_params['gi'].value * 1e-9 # to convert from nS to S
                    , 'max_conduction_speed': best_params['max_conduction_speed'].value, 
                    'poisson_rate': best_params['poisson_rate'].value,
                    'connector': best_params['connector'].value}
    netsim_obj.set_params(param_dict)
    print(f" best params {param_dict}")
    #create a dir for the output based on param space
    outputdir = os.path.expanduser('~') + '/temp/BEST/'
    netsim_obj.outputdir = outputdir
    #run the simulation
    r = netsim_obj.run()
    netsim_obj.get_results(reshape=True, reshape_spikes=True)
    matlab_params = netsim._add_matlab_params(netsim_obj.params)
    lfp = nslfp(matlab_params, netsim_obj.voltages, netsim_obj.ge, netsim_obj.gi)
    #plot the results
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=False)
    axes[0].plot(lfp[0])
    axes[0].set_xlabel('Time (ms)')
    axes[0].set_ylabel('LFP (V)')
    axes[1].eventplot(netsim_obj.spiketimes)
    axes[1].set_xlabel('Time (ms)')
    axes[1].set_ylabel('Neuron ID')
    f, Pxx = signal.welch(lfp[0], fs=1/matlab_params['dt'], nperseg=512, nfft=2048)
    axes[2].plot(f, Pxx)
    axes[2].set_xlim([0, 500])
    axes[2].set_xlabel('Frequency (Hz)')
    axes[2].set_ylabel('Power (V^2/Hz)')
    fig.savefig('best.png')
    plt.show()


if __name__ == '__main__':
    main()