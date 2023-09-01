# set environment for LAMMPS and msi2lmp executables
# to find potential and force field files
if ( "$?LAMMPS_POTENTIALS" == 0 ) setenv LAMMPS_POTENTIALS /home/kc32/.local/share/lammps/potentials
if ( "$?MSI2LMP_LIBRARY" == 0 ) setenv MSI2LMP_LIBRARY /home/kc32/.local/share/lammps/frc_files
