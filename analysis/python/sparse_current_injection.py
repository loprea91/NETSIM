import numpy as np
import argparse

#Only showing changes in values for neurons to keep sparse representation
#Assumes 1-1 mapping of stimulus to neuron

#file format is:
# uint32   - number of entries
# number of entry sets of the following one after another
# uint32   - start timestep
# uint32   - end timestep
# uint32   - number of neuron ids
# uint32[] - neuron ids
# float64  - stimulation intensity

parser = argparse.ArgumentParser(description='')
parser.add_argument('--dt', help='Timestep', required=True)
parser.add_argument('--output', help='Output matrix file', required=True)
args = parser.parse_args()

dt = float(args.dt)

#start step, end step, range of neurons to inject, intensity of injection
#assume that for a given timestep/intensity pair all stimulated neurons are included
values = [(0.2, 0.8, range(0,     50000),   5e-9),
          (0.0, 1.0, range(50000, 100000),  10e-9)
          ]

f = open(args.output, "wb")

#write number of timestep/neuron changes
f.write(np.uint32(len(values)).tobytes(order="C"))
for v in values:
    #start step value
    f.write(np.uint32(v[0]/dt).tobytes(order="C"))

    #end step value
    f.write(np.uint32(v[1]/dt).tobytes(order="C"))
    
    #number of indices
    f.write(np.uint32(len(v[2])).tobytes(order="C"))
    
    #neuron indices
    line = np.array(v[2], dtype=np.uint32, order="C")
    f.write(line.tobytes(order="C"))

    #stimulation intensity
    f.write(np.float64(v[3]).tobytes(order="C"))
    
f.close()
