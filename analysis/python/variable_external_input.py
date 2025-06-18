import numpy as np
import argparse
import array as arr

parser = argparse.ArgumentParser(description='')
parser.add_argument('-T', help='End time', required=True)
parser.add_argument('-dt', help='Timestep size', required=True)
parser.add_argument('-N', help='Number of Neurons', required=True)
parser.add_argument('--output', help='Output matrix file', required=True)
args = parser.parse_args()

N = int(args.N)
timesteps = int(float(args.T)/float(args.dt))

matrix = np.zeros((N, timesteps), dtype=np.float64, order="C")

matrix[:,200:400] = 9;

f = open(args.output, "wb")
f.write(matrix.tobytes(order="C"))
f.close()
