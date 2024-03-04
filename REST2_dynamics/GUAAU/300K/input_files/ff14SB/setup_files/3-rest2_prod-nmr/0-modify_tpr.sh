#!/bin/bash

source ~/Soft/Gromacs/gmx_2018.8_gtx/bin/GMXRC

RUN="nmr" # "nmr" or "simrna"

#read num
nrep=$1

for((i=0;i<nrep;i++))
do
  # include contraint definition between atoms 80 and 82.
 sed -e '/; distance restraints for WC_stem.*/i ; Constrain bond between atoms 80 and 82\n#ifdef CONSTR\n#include "c6_c5.itp"\n#endif\n' ../1-scaled_topols-${RUN}/topol-nCtoCn-${i}.top > topol-constr-${i}.top
	
done
