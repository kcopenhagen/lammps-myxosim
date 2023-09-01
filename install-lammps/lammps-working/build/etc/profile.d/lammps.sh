# set environment for LAMMPS and msi2lmp executables
# to find potential and force field files
LAMMPS_POTENTIALS=${LAMMPS_POTENTIALS-/home/kc32/.local/share/lammps/potentials}
MSI2LMP_LIBRARY=${MSI2LMP_LIBRARY-/home/kc32/.local/share/lammps/frc_files}
export LAMMPS_POTENTIALS MSI2LMP_LIBRARY
