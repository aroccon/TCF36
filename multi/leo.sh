module load profile/candidate
module load nvhpc/25.3
module load hpcx-mpi/2.19
cp Makefile_leonardo Makefile
make clean
make
#mpirun -np 2 ./mhit36 
