#!/bin/bash
DIR=$(pwd)
mkdir -p $HOME/software
cp -r lammps-working $HOME/software/
cp remake.sh $HOME/software/
cp lammps_remake.sh $HOME/software/

cd $HOME/software/lammps-working

./remake.sh
cd $DIR
