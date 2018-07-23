#!/bin/bash

nargs=3

print_usage()                                                                                  
{
   echo "Usage: $0 <X> <Y> <Z>"
   echo " <X>: geometric model"
   echo " <Y>: attribute file"
   echo " <Z>: number of processors"
   echo " Note: simmodsuite-9.0-140906 should be specified in your path if not already done"
   echo " Note: openmpi should also be specified if not already done"
}  

if [ "$#" -ne "$nargs" ] ; then
	print_usage
        exit 1
fi

modelfile=$1
attribfile=$2
nprocs=$3

echo
echo "Running parallel BL mesher (make sure to check BLMesher.log)..."

if [[ $modelfile == *.x_t ]] || [[ $modelfile == *.xmt_txt ]] || [[ $modelfile == *.X_T ]] || [[ $modelfile == *XMT_TXT ]] ; then
	executable=BLMesher-openmpi-O-parasolid-9.0-140906
elif [[ $modelfile == *.sat ]] || [[ $modelfile == *.SAT ]] ; then
	executable=missing
elif [[ $modelfile == *.smd ]] || [[ $modelfile == *.SMD ]] ; then
	executable=BLMesher-openmpi-O-geomsim-9.0-140906
else
	echo File extension unknown for $modelfile
	exit 1
fi

mpirun -np $nprocs /Install/develop/Meshing/BLMesherParallel/bin/x86_64_linux/$executable $modelfile mesh.sms 1 0.00032 $attribfile > BLMesher.log 2>&1

