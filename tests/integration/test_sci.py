import sys

sys.path.append( './python_interface/' )

import os
import pyNetsim.netsim as netsim
import argparse
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

netsim_dir = "./"
param_file = netsim_dir + "tests/test15.parameters"
output_dir = "./data"
matrix_file = netsim_dir + "sci.dat"

#create netsim object
net = netsim.netsim(params=param_file, output_dir=output_dir, netsim_dir=netsim_dir,
                    compile_on_run=True, verbose=True)

net.set_param("$sparse_current_injection_path", matrix_file)

dt = net.params["dt"]

values = [(0.2, 0.8, range(0, 50000),  5e-9),
          (0.0, 1.0, range(50000, 100000),  10e-9)
          ]

f = open(matrix_file, "wb")

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

#run netsim
r = net.run()
if r != 0:
    print("Error with NETSIM")
    exit(-1)
results = net.get_results()
