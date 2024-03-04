#!/bin/bash

RUN="nmr" # "nmr" or "simrna"

#read num
nrep=$1

# creates "nrep" nvt run files from the "nrep" topol$i.top files
for((i=0;i<nrep;i++))
do
	gmx grompp -o topol$i.tpr -f nvt-equil_rest2.mdp -p ../1-scaled_topols-${RUN}/topol-nCtoCn-$i.top -c ../${RUN}-seed.gro -n ../${RUN}-index.ndx -r ../${RUN}-seed.gro
done

touch plumed.dat
