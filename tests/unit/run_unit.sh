#!/bin/bash

# Compile the source code
echo "testing std_rng..."
gcc test_rng.c ../../src/rng_dist.c ../../src/isaac64.c ../../src/rng.c -I ../../src/ -lm -o test_rng

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
    echo "Running test_rng..."
    ./test_rng
else
    echo "Compilation failed"
fi

#try with fast_binomial=True
gcc test_rng.c ../../src/rng_dist.c ../../src/isaac64.c ../../src/rng.c -I ../../src/ -lm -DFAST_BINOMIAL=True -o test_rng

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
    echo "Running test_rng..."
    ./test_rng
else
    echo "Compilation failed"
fi