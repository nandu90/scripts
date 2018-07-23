#!/bin/bash

# Specify the number of proc
nproc=16

# Debug
#mpirun -np $nproc /Install/develop-scorec/build-core-sim-gnu-dbg/test/chef  2>&1 | tee chef.log

# Opt
mpirun -np $nproc /home/jfang3/develop-scorec/build-core-sim-gnu-opt/test/chef  2>&1 | tee chef.log


