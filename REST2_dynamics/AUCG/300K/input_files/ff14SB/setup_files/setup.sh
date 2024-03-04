#!\bin\bash

# RUN="nmr"  # "nmr" or "simrna"
# 
# 1. to create scaled topologies
# cd 1-scaled_topols-${RUN}
# 
# vi 1-scale_topols_parmed.py  #change "ion_idxs" to match those of "Ks" atoms in "-seed.gro" file 
# 
# python 1-scale_topols_parmed.py > scaling_factors.txt
# 
# ./2-modify_topols.sh 20
# 
# 2. to setup/run REST2 equilibration run
# cd ../2-nvt_equil-${RUN}
# 
# ./1-NVT_generate_tpr.sh 20
# 
# sbatch 2-NVT_mpi.sh  #run on HPCs
# 
# 
# 3. to setup/run REST2 production run
# cd ../3-rest2_prod-${RUN}
# 
# ./0-modify_tpr.sh 20
# 
# ./1-REST_generate_tpr.sh 20
# 
# sbatch 2a-rest2-mpi.sh  #run on HPCs
# 
# sbatch 2b-rest2-mpi-extend.sh #run on HPCs

