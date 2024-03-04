#!/bin/bash

RUN="nmr" # "nmr" or "simrna"

title=$(head -n 1 ../${RUN}-seed.gro)
system1="RNA"

#echo Please give number of replicas
#read num
nrep=$1

for((i=0;i<nrep;i++))
do
# change title of the topology file
sed "s/Generic title/${title}/g" topol$i.top > topol_title$i.top
# change the name of system1 to a more appropriate one.
sed "s/system1/${system1}/g" topol_title$i.top > topol_mol$i.top
# change all occurrences of string 2C to C2 in the topol$i.top  
sed 's/2C/C2/g' topol_mol$i.top > topol_edit_$i.top
# change all occurrences of string 3C to C3 in	the topol$i.top
sed 's/3C/C3/g' topol_edit_$i.top > topol-nCtoCn-$i.top

done

rm topol_edit* topol_title* topol_mol*
