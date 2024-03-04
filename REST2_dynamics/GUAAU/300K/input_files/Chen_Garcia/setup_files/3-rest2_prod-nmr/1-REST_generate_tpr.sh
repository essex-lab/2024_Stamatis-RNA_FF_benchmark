#!/bin/bash

source ~/Soft/Gromacs/gmx_2018.8_mpi_plmd/bin/GMXRC

RUN="nmr" # "nmr" or "simrna"

#read num
nrep=$1

# creates nrep md run files from the nrep topologies
for((i=0;i<nrep;i++))
do
  # prepare tpr file
	gmx_mpi grompp -o topol$i.tpr -f md_rest2.mdp -p topol-constr-$i.top -c ../2-nvt_equil-${RUN}/confout$i.gro -t ../2-nvt_equil-${RUN}/state$i.cpt -n ../${RUN}-index.ndx -r ../${RUN}-seed.gro
done

touch plumed.dat
