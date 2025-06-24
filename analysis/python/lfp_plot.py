import argparse
from pathlib import Path
import array as arr
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser(description='Generate plot of the LFP at a time from NETSIM')
parser.add_argument('--data', help='folder with binary output files from NETSIM')
parser.add_argument('--time', help='time of network to plot')
args = parser.parse_args()

params_path = list(Path(args.data).glob("*params.txt"))[0]
params_f = open(params_path, "r")
for param_line in params_f.readlines():
    param_line = param_line.strip()
    if (param_line == "" or
        param_line[0] == "#"):
        continue
    param = param_line.split("=")[0].strip()
    value = param_line.split("=")[1].strip()
    if param == "N":
        N = int(value)
    elif param == "dt":
        dt = float(value)
    elif param == "Ne":
        Ne = int(value)
    elif param == "NrowE":
        NrowE = int(value)
    elif param == "Ee":
        Ee = float(value)
    elif param == "Ei":
        Ei = float(value)
    elif param == "T":
        T = float(value)
    elif param == "record_downsample_factor":
        downsample_factor = int(value)
        
time = int(float(args.time)/(dt * downsample_factor))

vm_path = list(Path(args.data).glob("*vm.bin"))[0]
vm = np.fromfile(vm_path, dtype=np.float32).reshape(-1, N)

ge_path = list(Path(args.data).glob("*ge.bin"))[0]
ge = np.fromfile(ge_path, dtype=np.float32).reshape(-1, N)

gi_path = list(Path(args.data).glob("*gi.bin"))[0]
gi = np.fromfile(gi_path, dtype=np.float32).reshape(-1, N)

vm = vm[time, 0:Ne]
ge = ge[time, 0:Ne]
gi = gi[time, 0:Ne]

tau = int(0.006/(dt * downsample_factor))
exc = np.roll(ge * (Ee - vm), tau)
inh = gi * (Ei - vm)
lz = stats.zscore(abs(exc - 1.65*inh)).reshape(NrowE, NrowE)

plt.figure(figsize=(6.4, 4.8), dpi=300)
plt.imshow(lz)
plt.title("NETSIM LFP Plot")
plt.xlabel("Neuron x-position")
plt.ylabel("Neuron y-position")
cbar = plt.colorbar(label="z-score")
plt.tight_layout()
plt.savefig("./lfp_plot.png")
plt.show()
