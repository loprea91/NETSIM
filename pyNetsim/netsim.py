import subprocess
import os
import numpy as np
import argparse
import math
import logging
from pathlib import Path
import hashlib
from . import utility

#create a logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#Find some folders that we will need
#this is the folder where this file is located, i.e. the python_interface folder
basefolder = os.path.dirname(os.path.abspath(__file__))

#netsim src folder should be up one level
netsimfolder = os.path.join(basefolder, '../')
if not os.path.exists(netsimfolder):
    logger.error('netsim src folder not found with default path: {}'.format(netsimfolder))
else:
    logger.debug('netsim src folder found with default path: {}'.format(netsimfolder))

class netsim():
    def __init__(self,
                 params=None,
                 compilekwargs=None,
                 runkwargs=[],
                 output_dir='./', 
                 netsim_dir='./',
                 compile_on_run=True,
                 verbose=True,
                 debug=False):
        """ A class to run the NETSIM model. This class is a wrapper around the NETSIM model, which is written in C. 
        handles the compilation and running of the model, as well as the parsing of the output files. 

        Args:
            params (dict or str): a dict of parameters to be passed to the model, or a path to a parameter file
            compilekwargs (list, optional): A list of args to be passed to the compiler/make call. Defaults to ['COMPILE_TYPE=performance'].
            runkwargs (list, optional): A list of args to be passed to netsim. Defaults to [].
            output_dir (str, optional): Output directory of result files. Defaults to '.'
            netsim_dir (std, optional): Directory where the root netsim src tree is.
            compile_on_run (bool, optional): Compile netsim once run is called, alternatively .compile() should be called seperately. Defaults to True.

        Raises:
            TypeError: if params is not a dict or str
        """
        self.params = {}
        if compilekwargs is None:
            self.compilekwargs = ['COMPILE_TYPE=performance']
        else:
            self.compilekwargs = compilekwargs
            
        self.runkwargs = runkwargs
        self.compile_on_run = compile_on_run
        self.last_run_status = -999
        # Create a temp directory for the output and netsim src
        self.outputdir = output_dir
        self.netsimdir = netsim_dir
        self.analyze_results = None  # for the analyze method
        
        if params is None:  # If the user doesn't specify params, load the example params
            # Load the default params
            self.set_params(_load_example_params())
            self.paramfile = self.tempparamfile
        elif isinstance(params, str):  # If the user specifies a path to a paramfile, load it
            self.paramfile = params
            self.set_params(parse_paramfile(params))
        elif isinstance(params, dict):  # If the user specifies a dict of params, set them
            self.set_params(params)
            self.paramfile = self.write_paramfile(params)
        else:
            raise TypeError('params must be a dict or a string')

        # Set verbose and configure stdout and stderr
        self.verbose = verbose
        if self.verbose:
            self.stdout = None
            self.stderr = None
        else:
            self.stdout = subprocess.DEVNULL
            self.stderr = subprocess.DEVNULL

        self.debug = debug
        self.spikes_ready = False
        self.continuous_ready = False
        self.position_ready = False

        if self.debug:
            logger.debug('netsim object created, with params: {}'.format(self.params))
            logger.debug('netsim dir specified as: {}'.format(self.netsimdir))
            logger.debug('output dir specified as: {}'.format(self.outputdir))

    def compile(self):
        #make sure flags are correct
        self.determine_flags_from_params()
        
        #compile the model
        #rather than changing working directory, we just call make with the -C flag
        subprocess.run(['make', '-C', self.netsimdir, 'clean'],
                       stdout=self.stdout, stderr=self.stderr, check=True)      
        subprocess.run(['make', *self.compilekwargs, '-C', self.netsimdir],
                       stdout=self.stdout, stderr=self.stderr, check=True)
            
    def run(self):
        if self.compile_on_run:
            self.compile()

        # Prepare log file in the output directory if verbose is False
        log_file_path = os.path.join(self.outputdir, "netsim_run.log") if not self.verbose else None
        log_file = None

        try:
            if not self.verbose:
                # Open the log file for writing if verbose is False
                log_file = open(log_file_path, "w")
                self.stdout = log_file
                self.stderr = log_file
            else:
                # Allow console output if verbose is True
                self.stdout = None
                self.stderr = None

            # Run the model
            r = subprocess.run([self.netsimdir + 'netsim',
                                '-f', self._write_temp_paramfile(),
                                '-o', self.outputdir,
                                *self.runkwargs],
                            stdout=self.stdout, stderr=self.stderr)
        finally:
            # Ensure the log file is closed if it was opened
            if log_file:
                log_file.close()

        self.hash_conn()
        self._clean_up()  # Remove the temp parameter file
        self.last_run_status = r
        return r
    
    def get_results(self, reshape=False, reshape_spikes=False, long_sim_file_id=None):
        #parse the binary output file(s) and return a dict of results
        if "output_file_code" not in self.params:
            self.params["output_file_code"] = 0

        file_code = str(round(float(self.params["output_file_code"]))).zfill(8)
        if long_sim_file_id != None:
            file_id = str(long_sim_file_id).zfill(4)
            spike_ids = os.path.join(self.outputdir, file_code + 'spk_ids_' + file_id + '.bin')
            spike_times = os.path.join(self.outputdir, file_code + 'spk_times_' + file_id + '.bin')
        else:
            spike_ids = os.path.join(self.outputdir, file_code + 'spk_ids.bin')
            spike_times = os.path.join(self.outputdir, file_code + 'spk_times.bin')
        
        #load the binary files
        self.spikeids = np.fromfile(spike_ids, dtype=np.int32)
        self.spiketimes = np.fromfile(spike_times, dtype=np.float64)
        self.spikes_ready = True
        if self.params['start_record_time'] == self.params['stop_record_time']: #no continuous results
            if reshape_spikes:
                self.spiketimes = self._reshape_spikes(self.spikeids, self.spiketimes)
            return {'spikeids': self.spikeids,
                    'spiketimes': self.spiketimes,
                    'voltages': None,
                    'ge': None,
                    'gi': None}
        
        #load the continuous results
        if long_sim_file_id != None:
            file_id = str(long_sim_file_id).zfill(4)
            voltages = os.path.join(self.outputdir, file_code + 'vm_' + file_id + '.bin')
            ge_file = os.path.join(self.outputdir, file_code + 'ge_' + file_id + '.bin')
            gi_file = os.path.join(self.outputdir, file_code + 'gi_' + file_id + '.bin')
        else:
            voltages = os.path.join(self.outputdir, file_code + 'vm.bin')
            ge_file = os.path.join(self.outputdir, file_code + 'ge.bin')
            gi_file = os.path.join(self.outputdir, file_code + 'gi.bin')
        
        self.voltages = np.fromfile(voltages, dtype=np.float32)
        self.ge = np.fromfile(ge_file, dtype=np.float32)
        self.gi = np.fromfile(gi_file, dtype=np.float32)
        self.continuous_ready = True
        #reshape the results
        if reshape:
            #only reshape continuous results
            self.voltages = self._reshape_results(self.voltages)
            self.ge = self._reshape_results(self.ge)
            self.gi = self._reshape_results(self.gi)
        if reshape_spikes:
            self.spiketimes = self._reshape_spikes(self.spikeids, self.spiketimes)

        #return the results as a dict
        return {'spikeids': self.spikeids,
                'spiketimes': self.spiketimes,
                'voltages': self.voltages,
                'ge': self.ge,
                'gi': self.gi}
    
    def analyze(self, function_names, dt=None, bin_size=None):
        """
        Analyze the simulation output using specified functions from utility module.

        Parameters:
            function_names (list): List of function names as strings to be executed, 
                                or a predefined suite like {'AI'}.
            dt (float, optional): Time step size for functions that require it.
            bin_size (int, optional): Bin size for functions that require it.

        Returns:
            dict: Dictionary with function names as keys and their results as values.
        """
        analysis_suites = {
            'AI': ['pairwise_spike_correlation', 'coefficients_of_variation', 'fano_factor', 'firing_rate', 'power_spectrum_analysis']
        }

        expanded_function_names = []
        for name in function_names:
            if name in analysis_suites:
                expanded_function_names.extend(analysis_suites[name])
            else:
                expanded_function_names.append(name)

        if not self.spikes_ready:
            self.get_results()

        # Create a raster only if needed
        if any(fn in ['sttc', 'power_spectrum_analysis', 'firing_rate_raster', 'fano_factor', 'coefficients_of_variation_raster', 'pairwise_spike_correlation_raster'] for fn in expanded_function_names):
            if not hasattr(self, 'raster'):
                N = int(self.params['N'])
                T = int(np.ceil(self.params['T'] / dt))
                self.raster = utility.make_raster({'spk_ids': self.spikeids, 'spk_times': self.spiketimes}, N, T, dt)

        results = {}
        for func_name in expanded_function_names:
            func = getattr(utility, func_name, None)
            if not func:
                print(f"Function {func_name} not found in utility module.")
                continue

            # Argument checks and function calls
            if func_name == 'sttc':
                if bin_size is None:
                    print("Bin size is required for STTC.")
                    continue
                result = func(self.raster, bin_size)
            elif func_name == 'power_spectrum_analysis':
                if dt is None:
                    print("Time step size (sampling frequency) is required for power spectrum analysis.")
                    continue
                result = func(self.raster, dt)
            elif func_name == 'firing_rate':
                result = func({'spk_ids': self.spikeids, 'spk_times': self.spiketimes}, self.params['T'])
            elif func_name == 'firing_rate_raster':
                if dt is None:
                    print("Time step size is required for firing rate calculation.")
                    continue
                result = func(self.raster, dt)
            elif func_name == 'fano_factor':
                if bin_size is None:
                    print("Bin size is required for Fano factor calculation.")
                    continue
                result = func(self.raster, bin_size)
            elif func_name == 'coefficients_of_variation_raster':
                if dt is None:
                    print("Time step size is required for coefficient of variation calculation.")
                    continue
                result = func(self.raster, dt)
            elif func_name == 'coefficients_of_variation':
                result = func({'spk_ids': self.spikeids, 'spk_times': self.spiketimes})
            elif func_name == 'pairwise_spike_correlation':
                if bin_size is None:
                    print("Bin size is required for pairwise spike correlation.")
                    continue
                if dt is None:
                    print("Time step size is required for firing rate calculation.")
                    continue
                result = func({'spk_ids': self.spikeids, 'spk_times': self.spiketimes}, bin_size, dt)
            elif func_name == 'pairwise_spike_correlation_raster':
                if bin_size is None:
                    print("Bin size is required for pairwise spike correlation.")
                    continue
                result = func(self.raster, bin_size)
            else:
                result = func(self.raster)

            results[func_name] = result

        self.analyze_results = results
        return results
    
    def metric(self, analyze_function_names, target_values, results=None, dt=None, bin_size=None):
        """
        Calculate the Euclidean distance between the mean of the analyze methods' outputs and the target values.

        Parameters:
        - analyze_function_names: List of function names from the analyze method to be executed.
        - target_values: List of float target values.
        - results: Optional dictionary of precomputed analyze results.
        - dt: Time step size for functions that require it (optional).
        - bin_size: Bin size for functions that require it (optional).

        Returns:
        - distance: Euclidean distance between the mean of the analyze outputs and the target values.
        """
        # Use provided results, stored results, or compute new results
        if results is not None:
            analyze_results = results
        elif self.analyze_results is not None:
            analyze_results = self.analyze_results
        else:
            analyze_results = self.analyze(analyze_function_names, dt=dt, bin_size=bin_size)

        if len(analyze_function_names) != len(target_values):
            raise ValueError("The number of function names must match the number of target values.")

        # Calculate the mean of each analyze function's output
        means = []
        for func_name in analyze_function_names:
            output = analyze_results.get(func_name)
            if output is None:
                raise ValueError(f"Output for function {func_name} not found.")
            means.append(np.nanmean(output))

        # Calculate the Euclidean distance
        distance = np.linalg.norm(np.array(means) - np.array(target_values))

        return distance
    
    def get_position(self):
        file_code = str(round(self.params["output_file_code"])).zfill(8)
        pos_path = os.path.join(self.outputdir, file_code + 'pos.bin')
        self.pos = np.fromfile(pos_path, dtype=np.float64)
        self.pos = np.reshape(self.pos, (int(self.params['N']), 3));
        self.position_ready = True;
    
    def hash_conn(self):
        if self.params["save_connectivity"] == 1:
            # save connectivity as a hash value for easy comparison
            file_code = str(round(self.params["output_file_code"])).zfill(8)
            send_file = os.path.join(self.outputdir, file_code + 'ii.bin')
            receive_file = os.path.join(self.outputdir, file_code + 'jj.bin')
            hash_file = os.path.join(self.outputdir, file_code + 'hash.txt')

            send = np.fromfile(send_file, dtype=np.uint32)
            receive = np.fromfile(receive_file, dtype=np.uint32)
            combined_bytes = np.concatenate((send, receive)).tobytes()
            hash_value = hashlib.sha256(combined_bytes).hexdigest()
        
            print(f"Network Fingerprint: {hash_value}")
            with open(hash_file, 'w') as hash_file:
                hash_file.write(hash_value)

    def export_to_mat(self, save_dir):
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        if self.spikes_ready != True:
            self.get_results()
        savemat(save_dir / "spike_ids.mat", {"spike_ids": self.spikeids})
        savemat(save_dir / "spike_times.mat",{"spike_times": self.spiketimes})

        if self.position_ready != True:
            self.get_position()
        savemat(save_dir / "pos.mat", {"pos": self.pos})

        if self.params['start_record_time'] != self.params['stop_record_time']:
            if self.spikes_ready != True:
                self.get_results()
            savemat(save_dir / "ge.mat", {"ge": self.ge})
            savemat(save_dir / "gi.mat", {"gi": self.gi})
            savemat(save_dir / "vm.mat", {"vm": self.voltages})
            
    def _reshape_results(self, data):
        """Internal function to reshape the results from the binary files into a 2D array
        Follows the convention of the matlab version of the model

        Args:
            data (numpy array): a 1D array of results loaded from binary file

        Returns:
            result (numpy array): a 2D array of results in the shape (num_neurons, num_timesteps)
        """
        if 'dim' not in self.params:
            data_de = [data[int(idx)::int(self.params['N'])] for idx in range(int(self.params['N']))]
            data = np.array(data_de)
        #legacy code
        elif self.params['dim'] == 1:
            data = np.reshape( data, (int(self.params['nbins']), -1) )
        elif self.params['dim'] == 2:
            data = np.reshape( data, (int(self.params['nbins_exc']), -1) )
        return data

    def _reshape_spikes(self, spike_ids, spikes):
        """Internal function to reshape the spikes, into a list of numpy arrays containing the spike times for each neuron
        Follows the convention of the matlab version of the model

        Args:
            spike_ids (numpy array): a 1D array of spike ids loaded from binary file
            spikes (numpy array): a 1D array of spike times loaded from binary file

        Returns:
            result (numpy array): a 2D array of results in the shape (num_neurons, spike_times)
        """
        data = []
        for i in np.unique(spike_ids):
            data.append(spikes[spike_ids == i])
        #turn into ragged array
        data = np.array(data, dtype=object)
        return data

    def set_params(self, params=None, **kwargs):
        """_summary_

        Args:
            params (_type_): _description_
            **kwargs: _description_
        """
        
        if params is not None:
            self.params.update(params)
        self.params.update(kwargs)
        #ALSO ADD IN some matlab params for compatibility?
        #self.params = _add_matlab_params(self.params)
        self._write_temp_paramfile(self.params)

    def set_param(self, key, value):
        """_summary_

        Args:
            param (_type_): _description_
            **kwargs: _description_
        """
        if (key and value) is not None:
            self.params[key] = value
            self._write_temp_paramfile(self.params)
        
    def get_params(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        return self.params

    def write_paramfile(self, params):
        """_summary_

        Args:
            params (_type_): _description_

        Returns:
            _type_: _description_
        """
        paramfile = os.path.join(self.outputdir, 'temp.parameters')
        write_paramfile(params, paramfile)
        return paramfile

    def determine_flags_from_params(self):
        """_summary_ 
        
        Returns:
            _type_: _description_
        """
        current_injection_params = ("current_injection_count", "$current_injection_path")
        current_injection_arg = "CURRENT_INJECTION=yes"
        sparse_current_injection_params = ("$sparse_current_injection_path")
        sparse_current_injection_arg = "SPARSE_CURRENT_INJECTION=yes"
        poisson_params = ("poisson_rate",)
        poisson_arg = "EXTERNAL_INPUT=yes"
        variable_poisson_params = ("$poisson_rate_path",)
        variable_poisson_arg = "VARIABLE_EXTERNAL_INPUT=yes"
        long_simulation_params = ("partial_recording_time",)
        long_simulation_arg = "LONG_SIMULATION=yes"

        flags = []
        for key in self.params:
            if key in current_injection_params:
                flags.append(current_injection_arg)
            if key in sparse_current_injection_params:
                flags.append(sparse_current_injection_arg)
            if key in poisson_params:
                flags.append(poisson_arg)
            if key in variable_poisson_params:
                flags.append(variable_poisson_arg)
            if key in long_simulation_params:
                flags.append(long_simulation_arg)
        self.compilekwargs += list(set(flags))
        
    
    def _load_and_parse_paramfile(self, paramfile):
        self.params = parse_paramfile(paramfile)
    
    def _write_temp_paramfile(self, params=None):
        if params is None:
            params = self.params
        self.tempparamfile = self.write_paramfile(params)
        return self.tempparamfile
    
    def _clean_up(self):
        try:
            os.remove(self.tempparamfile)
        except:
            pass
    
    def __del__(self):
        self._clean_up()

    #output dir tools
    def _create_output_dir(self):
        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)
    
    @property
    def outputdir(self):
        return self._outputdir
    
    @outputdir.setter
    def outputdir(self, outputdir):
        self._outputdir = outputdir
        self._create_output_dir()
    
    @outputdir.getter
    def outputdir(self):
        return self._outputdir

def main():
    parser = argparse.ArgumentParser(prog="pynetsim",
                                     description="Compile, run and output graphs from NETSIM")
    parser.add_argument("-f", "--pfile", help="Path to parameter file", required=True)
    parser.add_argument("-p", "--path", help="Path to root netsim directory", required=True)
    parser.add_argument("-o", "--output", help="Path to desired output directory", default=".")
    parser.add_argument("-a", "--nargs", help="Arguments to pass to NETSIM", nargs='*')
    parser.add_argument("-c", "--compile", help="Compile before running", action="store_true", default=True)
    parser.add_argument("-v", "--verbose", action="store_true", help="Print extra information")
    args = parser.parse_args()
        
    #hop to the netsim dir, compile, run, and get results
    net = netsim(params=args.pfile,
                 netsim_dir=args.path,
                 output_dir=args.output,
                 compile_on_run=args.compile,
                 verbose=args.verbose)
    net.run()
    net.get_results(reshape=True)
    import matplotlib.pyplot as plt
    plt.scatter(net.spiketimes, net.spikeids, color='k', s=1, alpha=0.5)
    plt.figure()
    plt.plot(net.ge[4500])
    plt.show()


####  utility functions
def parse_paramfile(paramfile):
    #parse a parameter file
    with open(paramfile, 'r') as f:
        lines = f.readlines()
    params = {}
    for line in lines:
        if line[0] == '#' or line[0] == '\n':
            continue
        else:
            key, val = line.split('=')
            try:
                val = float(val.strip())
            except ValueError:
                val = val.strip()
            params[key.strip()] = val
    return params

def write_paramfile(params, paramfile='netsim.parameters'):
    #write a parameter file
    with open(paramfile, 'w') as f:
        for key, val in params.items():
            f.write('{} = {}\n'.format(key, val))
    return paramfile

def _load_example_params():
    #load the example parameters
    paramfile = os.path.join(netsimfolder, 'tests/test1.parameters')
    params = parse_paramfile(paramfile)
    return params

def _add_matlab_params(params):
    """
    Parse the parameters from the parameter dictionary; 
    adds additional parameters as added by importparams.m
    takes
        netsim: netsim object
    returns
        params: dictionary of parameters
    """
    
    #for ease of use, try casting the params to the correct type
    for key in params:
        try:
            params[key] = float(params[key])
        except:
            pass

    if params['connector'] == 4:
        params['dim'] = 2
    else:
        params['dim'] = 1

    #bin_size in latest refactor is ignored so we need to set it to 1
    params['bin_size'] = 1
    
    params['Ne'] = math.floor(0.8*params['N'])
    params['Ni'] = params['N'] - params['Ne']

    params['bin_size_x'] = params['bin_size'] ** (1/params['dim']) # number of neurons in one dimension of one bin (if 1D network, same as bin_size)
    params['nbins'] = params['N'] / (params['bin_size'] ** params['dim']) # number of bins (exc+inh)
    params['nbins_exc'] = 0.8 * params['nbins'] # number of bins (exc)
    params['nbins_inh'] = params['nbins'] - params['nbins_exc'] # number of bins (inh)

    DT = params['dt']*params['record_downsample_factor']
    
    params['t'] = np.arange(0, params['T'], DT) # time vector
    return params

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
