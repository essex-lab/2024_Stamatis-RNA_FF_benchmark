# 2024_Stamatis-RNA_FF_benchmark

This repository contains the published data, Jupyter Notebooks and input data to generate the plots found in *"Benchmarking RNA all-atom force fields using a set of diverse hairpin motifs"* by Dimitrios Stamatis, Chandra S. Verma, and Jonathan W. Essex.

Working paper version: https://doi.org/10.26434/chemrxiv-2023-0j7p4

Github repo: https://github.com/essex-lab/2024_Stamatis-RNA_FF_benchmark

Full dataset: https://doi.org/10.5281/zenodo.10715633

## Directory contents:

`gmx_FFs/` contains the Gromacs-compatible force field libraries used in this study.

`REST2_dynamics/` contains simulation input files (including .trp files for each replica), output trajectories, and reference structures used during analysis (base-replica trajectories are omitted; full dataset can be found on Zenodo).

`Analysis/notebooks/` contains Jupyter notebooks to calculate eRMSDs, G-vectors, η/θ pseudotorsions, and specific internucleotide distances and torsions, as well as to perform cluster-based anaylsis.

`Analysis/{UUCG,AUCG,CUUGU,GUAAU}/Gvec_conformations` contain input to reproduce published analysis data (pre-calculated G-vectors are omitted; full dataset can be found on Zenodo).

## Dependencies:
Python packages:
- Numpy
- Pandas
- Matplotlib
- Scikit-learn
- MDTraj<br>

Other:
- Gromacs (2018.8 or later)
- ParmEd
- baRNAba<br>
