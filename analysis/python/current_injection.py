import numpy as np
import argparse
import array as arr

parser = argparse.ArgumentParser(description='')
parser.add_argument('-T', help='End time', required=True)
parser.add_argument('-dt', help='Timestep size', required=True)
parser.add_argument('--count', help='Number of stimuli', required=True)
parser.add_argument('--output', help='Output matrix file', required=True)
args = parser.parse_args()

count = int(args.count)
timesteps = int(float(args.T)/float(args.dt))

matrix = np.zeros((count, timesteps), dtype=np.float64, order="C")
matrix[:,0:10] = 1.0

f = open(args.output, "wb")
f.write(matrix.tobytes(order="C"))
f.close()
