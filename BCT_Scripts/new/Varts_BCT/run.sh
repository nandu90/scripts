#!/bin/bash

echo "Check merge.inp and vbct.inp files properly"

rm -rf bctFiles
rm -rf mergedFiles
rm -rf ../bctInput
./make1
./make2
./make3

./run1
./run2
./run3
