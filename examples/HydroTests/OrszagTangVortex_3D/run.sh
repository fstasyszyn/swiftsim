#!/bin/bash

 # Generate the initial conditions if they are not present.
if [ ! -e OrszagTangVortex.hdf5 ]
then
    echo "Generating initial conditions for the OrszagTang " \
         "example..."
    python makeIC.py
fi

# Run SWIFT
../../swift --hydro --threads=4 OrszagTangVortex.yml 2>&1 | tee output.log

