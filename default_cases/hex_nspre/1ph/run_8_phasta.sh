rm nohup.out
nohup /Install/MPICH2/bin/mpirun -np 1 /home/nsaini3/develop_orig/develop/phasta/phSolver/200_memLS/phSolver/bin/x86_64_linux-icc/phastaIC.exe-mpich2-O < /dev/null &
#tail -f nohup.out

