#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e glassCube_64.hdf5 ]
then
    echo "Fetching initial glass file for the BrioWu example..."
    ./getGlass.sh
fi
if [ ! -e BrioWu.hdf5 ]
then
    echo "Generating initial conditions for the Brio Wu shock example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 BrioWu.yml 2>&1 | tee output.log

#python plotSolution.py 1
