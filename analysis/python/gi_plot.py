import argparse
from pathlib import Path
import array as arr
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Generate plot of conductances from NETSIM')
parser.add_argument('--data', help='folder with binary output files from NETSIM')
parser.add_argument('--id', help='id of neuron to plot')
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

id = int(args.id)

gi_path = list(Path(args.data).glob("*gi.bin"))[0]
gi_f = open(gi_path, "rb")
gi = arr.array("f", gi_f.read())
gi_f.close()

time = int(len(gi)/N)
times = np.linspace(0, time * dt, time)

plt.figure(figsize=(6.4, 4.8), dpi=300)
plt.plot(times, gi[id*time:(id+1)*time])
plt.title("NETSIM Inhibitory Conductance Plot")
plt.xlabel("Time(s)")
plt.ylabel("Conductance(S)")
plt.tight_layout()
plt.savefig("./gi_plot.png")
plt.show()
