import argparse
from pathlib import Path
import array as arr
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Generate scatter plot of spikes from NETSIM')
parser.add_argument('--data', help='folder with binary output files from NETSIM', required=True)
args = parser.parse_args()

ids_path = list(Path(args.data).glob("*spk_ids.bin"))[0]
times_path = list(Path(args.data).glob("*spk_times.bin"))[0]

ids_f = open(ids_path, "rb")
times_f = open(times_path, "rb")

ids = arr.array("I", ids_f.read())
times = arr.array("d", times_f.read())

ids_f.close()
times_f.close()

plt.figure(figsize=(6.4, 4.8), dpi=300)
plt.scatter(times, ids, s=0.025, c="Blue", facecolors="Blue", alpha=0.5)
plt.title("NETSIM Spike Plot")
plt.xlabel("Time(s)")
plt.ylabel("Neuron ID")
plt.tight_layout()
plt.savefig("./spike_plot.png")
plt.show()
