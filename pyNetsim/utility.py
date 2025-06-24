import numpy as np
import os
from scipy.signal import welch
from scipy import sparse
import multiprocessing as mp
import warnings
import os
import numpy as np

def read_spike_data(base_directory: str):
    """
    Reads spike times and spike ids from the NETSIM output folders within the given directory.

    Parameters:
    base_directory (str): The base directory

    Returns:
    dict: Dictionary with fields:
        spk_ids: neuron id (zero-based indexing)
        spk_times: spike times (s)

    Lawrence Oprea 2024
    """
    spike_data = {}

    if not os.path.exists(base_directory):
        print(f"Folder {base_directory} does not exist.")
        return spike_data

    spk_ids_file = None
    spk_times_file = None

    # Walk through the directory tree to find the files
    for root, _, files in os.walk(base_directory):
        for file in files:
            if file.endswith("spk_ids.bin"):
                spk_ids_file = os.path.join(root, file)
            elif file.endswith("spk_times.bin"):
                spk_times_file = os.path.join(root, file)

        # If both files are found, break the loop
        if spk_ids_file and spk_times_file:
            break

    if not spk_ids_file:
        print(f"spk_ids.bin not found in {base_directory}")
        return spike_data

    if not spk_times_file:
        print(f"spk_times.bin not found in {base_directory}")
        return spike_data

    try:
        spk_ids = np.fromfile(spk_ids_file, dtype=np.uint32)
        spk_times = np.fromfile(spk_times_file, dtype=np.double)
        spike_data = {"spk_ids": spk_ids, "spk_times": spk_times}

    except Exception as e:
        print(f"Error reading spike data from {base_directory}: {e}")

    return spike_data

def make_raster(data, N, T, dt):
    """
    Create a binary raster matrix from lists of spike times and spike ids.
    
    Parameters:
    - data: Dictionary {"spk_ids": [int], "spk_times": [int]}
    - N: Total number of neurons (rows in the raster)
    - T: Number of time points (columns in the raster)
    - dt: Time step (to convert times to integers)
    
    Returns:
    - rasters: 2D numpy array (N x T) binary raster matrixes

    Lawrence Oprea 2024
    """

    # Initialize then set spikes
    raster = np.zeros((N, T), dtype=int)
    raster[data['spk_ids'], (data['spk_times'] * (1/dt)).astype(int)] = 1

    return raster

def bin_spike_train(raster, bin_size):
    """
    Bin the spike train data using a specified bin size.

    Parameters:
    - raster: 2D numpy array (N x T) where each row represents the spike train of a neuron
    - bin_size: Size of the time window over which to compute spike counts

    Returns:
    - binned_raster: 2D numpy array where each row contains the binned spike counts for a neuron

    Lawrence Oprea 2024
    """
    num_neurons, num_timepoints = raster.shape
    num_bins = num_timepoints // bin_size
    binned_raster = raster[:, :num_bins * bin_size].reshape(num_neurons, num_bins, bin_size).sum(axis=2)
    return binned_raster

def pairwise_spike_correlation(data, bin_size, dt):
    """
    Compute the pairwise spike correlation efficiently using sparse matrices.
    
    Parameters:
      - data: Dictionary {"spk_ids": [int], "spk_times": [float]}.
      - bin_size: Number of time steps per bin.
      - dt: Time step (s); used to compute the bin duration.
    
    Returns:
      - neuron_ids: Array of unique neuron IDs included in the correlation matrix.
      - correlation_matrix: Pairwise correlation matrix.
      - decimation_factor: The decimation factor applied (default is 10, increased if memory errors occur).
      
    Note:
      This function bins the spike times into a spike count matrix, excludes neurons with zero variance,
      and then computes the Pearson correlation between the binned spike counts. If a MemoryError 
      occurs during correlation computation, the function decimates the spike count matrix by dropping rows.
    
    Lawrence Oprea 2024
    """

    import numpy as np
    from scipy.sparse import coo_matrix

    spk_ids = data['spk_ids']
    spk_times = data['spk_times']

    # Map neuron IDs to indices
    neuron_ids, neuron_id_indices = np.unique(spk_ids, return_inverse=True)

    # Define bin edges based on bin size and dt
    bin_duration = bin_size * dt
    min_time = spk_times.min()
    max_time = spk_times.max()
    num_bins = int(np.ceil((max_time - min_time) / bin_duration))
    bin_edges = min_time + np.arange(num_bins + 1) * bin_duration

    # Assign spikes to bins
    bin_indices = np.searchsorted(bin_edges, spk_times, side='right') - 1
    bin_indices = np.clip(bin_indices, 0, num_bins - 1)

    # Create the spike count matrix
    spike_count_matrix = coo_matrix(
        (np.ones_like(spk_times), (neuron_id_indices, bin_indices)),
        shape=(len(neuron_ids), num_bins)
    ).toarray()

    # Exclude neurons with zero variance
    variances = spike_count_matrix.var(axis=1)
    valid_indices = variances > 0
    spike_count_matrix = spike_count_matrix[valid_indices]
    neuron_ids = neuron_ids[valid_indices]

    # Compute correlation matrix
    decimation_factor = 10
    while True:
        try:
            correlation_matrix = np.corrcoef(spike_count_matrix)
            return neuron_ids, correlation_matrix, decimation_factor
        except MemoryError:
            decimation_factor *= 2
            if spike_count_matrix.shape[0] // decimation_factor < 2:
                return neuron_ids, np.array([[]]), decimation_factor
            spike_count_matrix = spike_count_matrix[::decimation_factor, :]


def pairwise_spike_correlation_raster(raster, bin_size):
    """
    Calculate pairwise spike correlation from the raster plot.
    
    Parameters:
    - raster: 2D numpy array (N x T) binary raster matrix
    - bin_size: Size of the bin for correlation calculation
    
    Returns:
    - correlations: 2D numpy array (N x N) of pairwise correlations
      Note: output matrix is symmetric and likely contains many NaNs

    Lawrence Oprea 2024
    """

    binned_raster = bin_spike_train(raster, bin_size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        correlations = np.corrcoef(binned_raster)
    
    return correlations

def coefficients_of_variation_raster(raster, timestep):
    """
    Calculate the coefficient of variation of interspike intervals for each neuron.
    
    Parameters:
    - raster: 2D numpy array (N x T) binary raster matrix
    - timestep: Time step size
    
    Returns:
    - cv: 1D numpy array of coefficients of variation for each neuron

    Lawrence Oprea 2024
    """
    num_neurons = raster.shape[0]
    cv = np.full(num_neurons, np.nan)  # Initialize with NaN
    
    for i in range(num_neurons):
        spike_times = np.where(raster[i] == 1)[0] * timestep
        if len(spike_times) > 1:
            isi = np.diff(spike_times)
            mean_isi = np.mean(isi)
            if mean_isi != 0:
                cv[i] = np.std(isi) / mean_isi
    
    return cv

import numpy as np

def coefficients_of_variation(data):
    """
    Calculate the coefficient of variation of interspike intervals (ISIs) for each neuron.

    Parameters:
        data (dict): A dictionary containing:
            'spk_ids' (array): An array of neuron IDs (0-indexed).
            'spk_times' (array): An array of corresponding spike times (in seconds).

    Returns:
        cv (np.array): A NumPy array of the same length as the number of unique neuron IDs 
                       in 'spk_ids'. Each element is the coefficient of variation for a neuron.
                       Neurons with fewer than 2 spikes will have NaN values, as will neurons
                       that have a zero mean ISI.
    
    Lawrence Oprea 2025
    """

    spk_ids = data['spk_ids']
    spk_times = data['spk_times']
    unique_ids = np.unique(spk_ids)
    
    cv = np.full(len(unique_ids), np.nan)

    for i, neuron_id in enumerate(unique_ids):
        spike_times_neuron =  np.sort(spk_times[spk_ids == neuron_id])

        if len(spike_times_neuron) >= 2:
            isi = np.diff(spike_times_neuron)
            mean_isi = np.mean(isi)
            if mean_isi != 0:
                cv[i] = np.std(isi) / mean_isi

    return cv

def firing_rate(data, T):
    """
    Calculate the firing rate in Hz for each neuron.
    
    Parameters:
    - data: dictionary of spike ids and times 
    - T: Total time (s)
    
    Returns:
    - 1D numpy array of firing rates for each neuron

    Lawrence Oprea 2025
    """
    spk_ids = data["spk_ids"]
    spike_counts = np.bincount(spk_ids)
    firing_rates = spike_counts / T

    return firing_rates

def firing_rate_raster(raster, dt):
    """
    Calculate the firing rate in Hz for each neuron.
    
    Parameters:
    - raster: 2D numpy array (N x T) binary raster matrix
    - dt: Time step size
    
    Returns:
    - 1D numpy array of firing rates for each neuron

    Lawrence Oprea 2024
    """
    
    return np.mean(raster, axis=1) / dt

def fano_factor(raster, bin_size):
    """
    Calculate the Fano factor for each neuron's spike train using a specified bin size.

    Parameters:
    - raster: 2D numpy array (N x T) where each row represents the spike train of a neuron
    - bin_size: Size of the time window over which to compute spike counts

    Returns:
    - fano_factors: 1D numpy array (N,) containing the Fano factor for each neuron
      Note: Fano factors are NaN where mean spike count is zero
    
    Lawrence Oprea 2024
    """
    binned_raster = bin_spike_train(raster, bin_size)
    mean_counts = binned_raster.mean(axis=1)
    variance_counts = binned_raster.var(axis=1)
    fano_factors = np.where(mean_counts > 0, variance_counts / mean_counts, np.nan)
    return fano_factors

def power_spectrum_analysis(raster, fs):
    """
    Perform a power spectrum analysis on spike trains using Welch's method.

    Parameters:
    - raster: 2D numpy array (N x T)
        Binary spike matrix where each row represents the spike train of a neuron.
    - fs: float
        Sampling frequency of the spike train data (in Hz).

    Returns:
    - frequencies: 1D numpy array
        Array of frequency bins.
    - power_spectra: 2D numpy array (N x F)
        Power spectra for each neuron, where each row corresponds to a neuron's power spectral density.
    """
    # Compute power spectra for all neurons simultaneously
    frequencies, power_spectra = welch(raster, fs=1.0/fs, axis=1, nperseg=min(256, raster.shape[1]))
    return frequencies, power_spectra

def compute_P(spike_times_1, spike_times_2, bin_size):
    """
    Helper function for sttc.

    Lawrence Oprea 2024
    """
    if spike_times_1.size == 0 or spike_times_2.size == 0:
        return 0.0

    indices = np.searchsorted(spike_times_2, spike_times_1)
    count = 0
    for idx, t in zip(indices, spike_times_1):
        # Check spikes within bin_size on both sides
        window_indices = spike_times_2[max(0, idx - bin_size): min(len(spike_times_2), idx + bin_size + 1)]
        if np.any(np.abs(window_indices - t) <= bin_size):
            count += 1
    return count / len(spike_times_1)

def compute_sttc_pair(args):
    """
    Helper function for sttc.

    Lawrence Oprea 2024
    """
    i, j, spike_times, TA, bin_size = args
    num_spikes_i = len(spike_times[i])
    num_spikes_j = len(spike_times[j])

    if num_spikes_i == 0 or num_spikes_j == 0:
        return (i, j, 0.0)

    PA = compute_P(spike_times[i], spike_times[j], bin_size)
    PB = compute_P(spike_times[j], spike_times[i], bin_size)

    TA_i = TA[i]
    TB_j = TA[j]

    denom_A = 1 - PA * TB_j
    denom_B = 1 - PB * TA_i

    term_A = (PA - TB_j) / denom_A if denom_A != 0 else 0.0
    term_B = (PB - TA_i) / denom_B if denom_B != 0 else 0.0

    sttc_value = 0.5 * (term_A + term_B)
    return (i, j, sttc_value)

def sttc(raster, bin_size, use_mp=False):
    """
    Compute the Spike Time Tiling Coefficient (STTC) for each pair of neurons using a specified time window.

    Parameters:
    - raster: 2D numpy array (N x T) or scipy.sparse matrix
        Binary spike matrix where each row represents the spike train of a neuron.
    - bin_size: int
        Time window (in number of bins) within which to consider spikes "synchronized".
    - use_mp: bool
        Whether to use multiprocessing.

    Returns:
    - sttc_matrix: 2D numpy array (N x N)
        Symmetric matrix of STTC values between neuron pairs.

    Requires:
    - compute_P()
    - compute_sttc_pair()

    Lawrence Oprea 2024
    """
    N, T = raster.shape

    # Convert to sparse matrix if not already
    if not sparse.issparse(raster):
        raster = sparse.csr_matrix(raster)

    # Precompute spike times and TA for each neuron
    spike_times = []
    TA = np.zeros(N)

    for i in range(N):
        # Get spike times for neuron i
        spike_times_i = raster[i].indices  # Non-zero indices in the sparse row
        spike_times.append(spike_times_i)

        if spike_times_i.size == 0:
            TA[i] = 0.0
            continue

        # Compute TA for neuron i
        time_within_bin = np.zeros(T, dtype=bool)
        for t in spike_times_i:
            start = max(0, t - bin_size)
            end = min(T, t + bin_size + 1)
            time_within_bin[start:end] = True
        TA[i] = np.sum(time_within_bin) / T

    neuron_pairs = [(i, j, spike_times, TA, bin_size) for i in range(N) for j in range(i + 1, N)]

    sttc_matrix = np.zeros((N, N))

    if use_mp:
        # Use multiprocessing
        with mp.Pool() as pool:
            for i, j, value in pool.imap_unordered(compute_sttc_pair, neuron_pairs):
                sttc_matrix[i, j] = value
                sttc_matrix[j, i] = value  # Symmetric matrix
    else:
        # Single-threaded computation
        for args in neuron_pairs:
            i, j, value = compute_sttc_pair(args)
            sttc_matrix[i, j] = value
            sttc_matrix[j, i] = value  # Symmetric matrix

    return sttc_matrix