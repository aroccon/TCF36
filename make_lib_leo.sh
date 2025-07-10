git clone https://github.com/NVIDIA/cuDecomp
cd cuDecomp
mkdir build
cd build
module load nvhpc/24.3
module load openmpi/4.1.6--nvhpc--24.3 
cmake ..
#enable nvshmem
#cmake -DCUDECOMP_BUILD_EXTRAS=1 -DCUDECOMP_ENABLE_NVSHMEM=1 ..
make -j
