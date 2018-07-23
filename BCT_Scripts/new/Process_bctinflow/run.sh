#!/bin/bash


#Directions for running the script:
#Have all bct.dat.* files in the folder where you run the script
#Have Coin.f and gen_xyzts.f in the folder
#Have compile_command and compile_xyzts in the folder
#That's all

echo "Generating myrank.txt and xyzts.dat"

mkdir bct
mv bct.dat.* bct/

ls bct > myrank.txt

sed -i "s/bct.dat.//" myrank.txt

./compile_command
./Runme_Coin

./compile_xyzts
./xyzts


