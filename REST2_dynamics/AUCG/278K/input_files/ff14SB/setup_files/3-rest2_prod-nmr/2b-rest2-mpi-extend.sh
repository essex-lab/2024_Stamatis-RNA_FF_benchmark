#!/bin/bash
#
#SBATCH --partition batch      # Run on the batch partition ('gpu' and 'gtx1080' are alternatives)
#SBATCH --job-name=2K-14SB
#SBATCH --ntasks-per-node=40     # Tasks per node
#SBATCH --nodes=10               # Number of nodes requested
#SBATCH --time=54:00:00

module load openmpi/3.0.0/gcc
module load gcc/6.4.0
module load cmake/3.10.2
module load fftw/3.3.7/gcc

# source the GROMACS environment
source /home/ds2e19/Soft/Gromacs/gmx_2018.8_mpi_plmd/bin/GMXRC
export PLUMED_KERNEL=/home/ds2e19/Soft/Plumed/plmd_2.5.3/lib/libplumed_2.5.3Kernel.so

# check hard ware environment
NNODE=10
NPROCPERNODE=40
NPROC=$[$NNODE * $NPROCPERNODE]
#NTOMP=5
#NPERNODE=$[$NPROCPERNODE / $NTOMP]
#NTMPI=$[$NNODE * $NPERNODE]
#export OMP_NUM_THREADS=$NTOMP

##################################################################################

nrepl=20 #number of REST2 replicas

#Extend the simulation
mpirun -np $NPROC gmx_mpi mdrun -plumed plumed.dat -v -multi $nrepl -s topol -replex 2000 -hrex -dlb no -maxh 54 -cpi state -append &>> run.out
