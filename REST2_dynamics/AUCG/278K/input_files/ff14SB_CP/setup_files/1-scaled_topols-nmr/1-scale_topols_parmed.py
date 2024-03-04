import parmed as pmd
import os
import copy
import numpy as np


############## Define stuff here ##############
RUN="nmr" # "nmr" or "simrna"
topolpath = "../topol-hmr-%s.top"%RUN
gropath = "../%s-seed.gro"%RUN

# Use defined position and distance restraints.
use_posre = False
posrepath = "../%s-posre.itp"%RUN
use_disre = True
disrepath = "../disre_WC_stem.itp"

# Parameters to control the scaling range (1...tmin/tmax) and progression step.
tmin = 278
tmax = 927
nrepl = 20

hot_idx_i = 5  # 0-based index of first hot residue
hot_idx_f = 8  # 0-based index of last hot residue

scale_ions = True
ion_idxs = [6832, 6833, 6846, 6851] #0-based residue indexes of hot ions.

scale_dihedrals = True
# atom types for definition of the unscaled dihedrals
unscaled_dihed_atypes = {'H', 'HA', 'H4', 'H5', \
                'C', 'CA', 'CB', 'CP', 'CQ', 'CS', 'C4', 'C5', \
                'N*', 'NA', 'NB', 'NC', 'N2', \
		'O'}

os.system("rm -rf topol*")

############## Initiate stuff according to user input ##############

#atoms with any type of scaling have underscored atom types.
unscaled_dihed_atypes_mod = [n + '_' for n in unscaled_dihed_atypes] 

# set of hot residue indices (usually a continuous region).
wanted_indices = {i for i in range(hot_idx_i, hot_idx_f + 1)}
# gather all atom types that take part in unscaled dihedrals.
unwanted_types = set(unscaled_dihed_atypes)
for e in unscaled_dihed_atypes_mod:
    unwanted_types.add(e)

#read distance and position restraint files.
if use_disre:
    with open(disrepath, 'r') as disrefile:
        disre = disrefile.readlines()
if use_posre:
#    with open(posrepath, 'r') as posrefile:
#        posre = posrefile.readlines()
    posre = '; Include Position restraint file\n#ifdef POSRES\n#include "posre.itp"\n#endif\n'

# Read the topologies and save an unscaled local copy 
f = pmd.load_file(topolpath, xyz=gropath)
f.save("topol_unscaled.top")

# Dicts to bind scaled dihedrals to their original ones.
dihedral_dict1 = {}
dihedral_dict2 = {}

# Define highest scaling factor and geometric progression step.
scaling_factor = 1.0
scaling_step = 1 / np.exp(np.log(tmax / tmin) / (nrepl - 1))


############## Scale charges, LJ epsilon, and proper dihedrals ##############
for c in range(nrepl):
    f = pmd.load_file(topolpath, xyz=gropath)
    print("Scaling with factor: %s" % scaling_factor)
    
    #scale charge and epsilon of hot atoms (PLUMED recipe).
    for residue in f.residues[hot_idx_i:hot_idx_f + 1]:
        for atom in residue:
            atom.charge *= scaling_factor ** 0.5
            atom.atom_type = copy.deepcopy(atom.atom_type)
            atom.atom_type.name = atom.atom_type.name.strip("_") + "_"
            atom.atom_type.epsilon *= scaling_factor
            atom.type = atom.atom_type.name

    if scale_dihedrals:
        #scale hot region dihedrals (mostly PLUMED recipe).
        for dihedral in f.dihedrals:
            # get residue indices of the dihedral atoms.
            res_indices = [dihedral.atom1.residue.idx, dihedral.atom2.residue.idx, dihedral.atom3.residue.idx,
                           dihedral.atom4.residue.idx]
            # get types of the dihedral atoms
            res_types = [dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, dihedral.atom4.type]

            # check the extend wanted_indices and unwanted_types intersect with dihedral residue indeces and nucleobase atom types respectively
            intersect_idx = [j for j in res_indices if j in wanted_indices]
            intersect_type = [k for k in res_types if k in unwanted_types]
            # Following PLUMED recipe for dihedral scaling and my own intuition for excluding dihedrals.
            # if dihedral in hot region (i.e. all atoms in hot residues) and dihedral not in nucleobase (i.e. <=2 atom in unwanted types)
            if (len(intersect_idx) > 3 and len(intersect_type) <= 2):
                # modify and store newly-encountered dihedral types in dictionary
                if dihedral not in dihedral_dict1.keys():
                    new_dihedral_type = copy.deepcopy(dihedral.type)

                    # we need to treat ListType and Type dihedrals separately.
                    if "List" in str(type(new_dihedral_type)):
                        for dtype in new_dihedral_type:
                            dtype.phi_k *= scaling_factor
                    else:
                        new_dihedral_type.phi_k *= float(scaling_factor)

                    f.dihedral_types += [new_dihedral_type]
                    dihedral_dict1[dihedral.type] = f.dihedral_types[-1]

                # causes this list to "claim" all of the items it contains and subsequently indexes all of its items \
                # (i.e. writes the modified values back to the original dihedrals)
                f.dihedral_types.claim()
                # bind the dihedral back to its original one?
                dihedral.type = dihedral_dict1[dihedral.type]

            # if 1-3 atoms in hot region, use squared lambda scaling (second argument should always be true for nucleotides).
            elif (len(intersect_idx) >= 1 and len(intersect_idx) <= 3) and len(intersect_type) <= 2:
                # modify and store newly-encountered dihedral types in dictionary
                if dihedral not in dihedral_dict2.keys():
                    new_dihedral_type = copy.deepcopy(dihedral.type)

                    # we need to treat ListType and Type dihedrals separately.
                    if "List" in str(type(new_dihedral_type)):
                        for dtype in new_dihedral_type:
                            dtype.phi_k *= float(scaling_factor) ** 0.5
                    else:
                        new_dihedral_type.phi_k *= float(scaling_factor) ** 0.5

                    f.dihedral_types += [new_dihedral_type]
                    dihedral_dict2[dihedral.type] = f.dihedral_types[-1]

                # causes this list to "claim" all of the items it contains and subsequently indexes all of its items \
                # (i.e. writes the modified values back to the original dihedrals)
                f.dihedral_types.claim()
                # bind the dihedral back to its original one?
                dihedral.type = dihedral_dict2[dihedral.type]

            # #else no scaling
            else:
                pass

    if scale_ions:
        #scale ion epsilon and charge (creates new molecule, atom and residue types; PLUMED recipe).
        for i in ion_idxs:
            ion = f.residues[i]
            ion.name = ion.name + "s"
            for atom in ion:
                atom.charge *= scaling_factor ** 0.5
                atom.atom_type = copy.deepcopy(atom.atom_type)
                atom.atom_type.name = atom.atom_type.name.strip("_") + "s_"
                atom.atom_type.epsilon *= scaling_factor
                atom.type = atom.atom_type.name
                atom.name = atom.name + "s"
    
    fname = "topol%d.top" % c

    f.save(fname)

    if use_posre or use_disre:
        # Insert distance restraints to the end of first moleculetype.
        #(may be extended to include position restraints and/or other moleculetypes)
        contents = open(fname).readlines()
        new_contents = [[]]
        #parse the [ ... ] fields of topol#.top in new_contents list.
        for line in contents:
            if line[0] == "[":
                new_contents += [[]]
            new_contents[-1] += [line]
        #append the distance resraint file content after the first [ dihedrals ] instance.
        dihedral_index = [i for i, x in enumerate(new_contents) if "[ dihedrals ]" in x[0]][0]
        if use_disre: new_contents.insert(dihedral_index + 1, disre + ["\n"])
        if use_posre: new_contents.insert(dihedral_index + 1, [posre] + ["\n"])
        #flatten the list and (over)write topol#.top
        new_lines = sum(new_contents, [])
        with open(fname, "w") as file:
            for line in new_lines:
                file.write(line)

    #update scaling factor for next replica.
    scaling_factor *= scaling_step

