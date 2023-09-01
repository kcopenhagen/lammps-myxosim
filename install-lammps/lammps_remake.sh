#!/bin/bash

VERSION=working
cd lammps-${VERSION}
cd build

module purge
module load intel/19.1.1.217
module load intel-mpi/intel/2019.7

cmake3 -D CMAKE_INSTALL_PREFIX=$HOME/.local \
-D LAMMPS_MACHINE=della \
-D ENABLE_TESTING=no \
-D BUILD_MPI=yes \
-D BUILD_OMP=yes \
-D CMAKE_BUILD_TYPE=Release \
-D CMAKE_CXX_COMPILER=icpc \
-D CMAKE_CXX_FLAGS_RELEASE="-Ofast -xHost -qopenmp -restrict -DNDEBUG" \
-D PKG_MOLECULE=yes \
-D PKG_BROWNIAN=yes \
-D PKG_EXTRA-FIX=yes \
-D PKG_DIPOLE=yes \
-D PKG_RIGID=yes \
-D PKG_KSPACE=yes -D FFT=MKL -D FFT_SINGLE=yes \
-D PKG_INTEL=yes -D INTEL_ARCH=cpu -D INTEL_LRT_MODE=threads ../cmake

make -j 10
make install
