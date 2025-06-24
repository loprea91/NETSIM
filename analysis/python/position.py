import argparse
from pathlib import Path
import array as arr
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Generate scatter plot of spikes from NETSIM')
parser.add_argument('--data', help='folder with binary output files from NETSIM', required=True)
args = parser.parse_args()

pos_path = list(Path(args.data).glob("*pos.bin"))[0]
pos_f = open(pos_path, "rb")
pos = arr.array("d", pos_f.read())
pos_f.close()

#make N by 3 Matrix
pos = np.array(pos).reshape((-1, 3))

#extract coordinates
X = pos[:,0]
Y = pos[:,1]
Z = pos[:,2]

fig = plt.figure(figsize=(6.4, 4.8), dpi=300)
ax = fig.add_subplot(projection='3d')

ax.scatter(X, Y, Z, s=0.25, c="Blue")
plt.title("NETSIM Neuron Positions")
ax.set_xlabel("X-axis (m)")
ax.set_ylabel("Y-axis (m)")
ax.set_zlabel("Z-axis (m)")
plt.savefig("./position_plot.png")
