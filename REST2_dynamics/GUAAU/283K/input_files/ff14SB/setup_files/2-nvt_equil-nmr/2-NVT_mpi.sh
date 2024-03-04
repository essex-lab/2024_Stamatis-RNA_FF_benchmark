#!/bin/bash
#
#SBATCH --job-name=2KC-14SB
#SBATCH --ntasks-per-node=40     # Tasks per node
#SBATCH --nodes=2                # Number of nodes requested
#SBATCH --time=04:00:00

module load openmpi/3.0.0/gcc
module load gcc/6.4.0
module load cmake/3.10.2
module load fftw/3.3.7/gcc

# source the GROMACS environment
source /home/ds2e19/Soft/Gromacs/gmx_2018.8_mpi_plmd/bin/GMXRC
export PLUMED_KERNEL=/home/ds2e19/Soft/Plumed/plmd_2.5.3/lib/libplumed_2.5.3Kernel.so

# check hard ware environment
NNODE=2 	 #should be the same as "nodes"
NPROCPERNODE=40	 #should be the same as "ntasks-per-node"
NTOMP=4          #OpenMP threads per task
NPERNODE=$[$NPROCPERNODE / $NTOMP] #MPI tasks per node
NTMPI=$[$NNODE * $NPERNODE]	   #total MPI tasks
export OMP_NUM_THREADS=$NTOMP

##################################################################################

nrepl=20

# Run NVT equilibration before H-REMD (REST2) for nrepl replicas with exchange 
# attempts every 5000000*0.002fs = 10ns (basically no exchange will be performed:P)
# replex is set to a number greater than the nvt run time, so no exchanges occur in equilibration
mpirun -np $NTMPI -npernode $NPERNODE gmx_mpi mdrun -plumed plumed.dat -ntomp $NTOMP -v -multi $nrepl -s topol -replex 5000000 -hrex
