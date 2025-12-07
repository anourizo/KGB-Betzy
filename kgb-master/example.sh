#!/bin/bash

#SBATCH --account=NN12015K
#SBATCH --job-name=KGB
#SBATCH --time=05:00:00

#SBATCH --nodes=10 --ntasks=1200 --cpus-per-task=1

##SBATCH --nodes=4 --ntasks=16 --cpus-per-task=1

#SBATCH --output=slurm-test-%J.out
#SBATCH --mail-type=FAIL

ml load foss/2022a GSL/2.7-GCC-11.3.0 HDF5/1.12.2-gompi-2022a

export LD_LIBRARY_PATH=$HOME/software/cfitsio/lib:$HOME/software/Healpix_3.83/lib:$LD_LIBRARY_PATH
##COMMAND##
time srun -n 1200 ./KGBevolution -n 40 -m 30 -s settings.ini

