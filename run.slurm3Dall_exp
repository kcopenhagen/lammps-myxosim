#!/bin/bash
#SBATCH --job-name=myxo-sim       # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=9              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=75M       # memory per cpu-core (4G is default)
#SBATCH --time=23:00:00          # total run time limit (HH:MM:SS)
#SBATCH --constraint=cascade,skylake

module purge
module load matlab/R2020b
matlab -singleCompThread -nodisplay -nosplash -r Initial_all

module purge
module load intel/19.1.1.217
module load intel-mpi/intel/2019.7
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun $HOME/.local/bin/lmp_della -in in.spfr3Dall_exp

