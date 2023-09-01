#!/bin/bash

NAME=kc32

DIR=$(pwd)
mkdir $HOME/software
cp -r lammps-working $HOME/software/
cp remake.sh $HOME/software/
cp lammps_remake.sh $HOME/software/

cd $HOME/software/build
sed s/kc32/$NAME/ cmake_install.cmake > cmaketemp
mv cmaketemp cmake_install.cmake

cd ../
./remake.sh
cd $DIR
