import argparse
from pathlib import Path
import array as arr
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

parser = argparse.ArgumentParser(description='Generate plot of connection matrix from NETSIM')
parser.add_argument('--data', help='folder with binary output files from NETSIM', required=True)
parser.add_argument('--output', help='folder to output file', default='.')
args = parser.parse_args()

output = Path(args.output)

ii_path = list(Path(args.data).glob("*ii.bin"))[0]
ii_f = open(ii_path, "rb")
ii = arr.array("I", ii_f.read())
ii_f.close()

jj_path = list(Path(args.data).glob("*jj.bin"))[0]
jj_f = open(jj_path, "rb")
jj = arr.array("I", jj_f.read())
jj_f.close()

#G = nx.Graph(zip(ii, jj))

plt.figure(figsize=(20, 16))
plt.scatter(ii, jj, s=20, alpha=1)
plt.title("NETSIM Connection Matrix")
plt.xlabel("Sending index")
plt.ylabel("Recieving index")
plt.savefig(output / "connection_matrix_plot.png")
