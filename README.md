
# EasyTrajH5

Trajectory management for mdtraj H5 files

- EasyH5File wrapper for h5py
  - methods-based interface for dataset/attr management
  - string, json, file import/export methods
  - useful schema dictionary
- select_mask with new atom selection language for MD
  - proper set operations using 'not', 'diff', 'merge', 'intersect'
  - allows both amber and mdtraj style selections
  - proper 0-based indexing using flexible numbering list and ranges
  - user defined residue classes
- EasyTrajH5File
  - drop-in replacement for mdtraj.H5TrajecotryFile
  - subclasses EasyH5File with all above methods
  - allows simple override to h5pyd
  - efficient streaming of mdtraj.Trajectory from file using fancy indexing
    and atom masks with select_mask
- conveient transforms between mdtraj, parmed and openmm
- cli for interrogating h5 files


### Selection Language

Selects atom based on a selection string which is based on a combination
atom selection langauge.There are several different types of selection
modes strings, where the first word is often used a mode selector.

- keywords
    - If more than one keyword is specified, it is assumed they are joined with "or"
      operation (i.e. 'ligand protein' will return both ligand and protein atom indices).
    - accepts (in any order): 'ligand', 'protein', 'water', 'lipid', 'salt',
      'solvent', 'lipid', 'nucleic', 'resname', 'resid', 'atom'
    - 'ligand' will find the residue 'LIG', 'UNL', 'UNK', or
       whatever is in 'ligand' in the h5 'easytrajh5/data/select.yaml'
    - 'pocket' will find the closest 6 residues to the 'ligand' group.
    - 'near' will require a following resname, with an optional integer, e.g.:
        'near ATP'
        'near ATP 5'
    - 'resname' identifies a single residue type (usually a ligand):
        'resname LEU'
    - 'resi' for residue 0-indexed selections
        "resi 0 10-13" - selects atoms in the first and 11th, 12th and 13th residues
    - 'atom' atom 0-indexed selections
        "atom 0 55 43 101-105" - selects the first, 56th, 44th, 101 to 105th atom
- AMBER-style atom selection https://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax
    "amber :ALA,LYS" - selects all alanine and lysine residues
- MDTraj-style atom selection - https://mdtraj.org/1.9.4/atom_selection.html
    "mdtraj protein and water" - selects protein and water
- furthermore, selections can be combined with set operators ("not", "intersect", "merge", "diff"),
    "intersect {not {amber :ALA}} {protein}"
    "diff {protein} {not {amber :ALA}}"
    "not {resname LEU}"
    "merge {near BSM 8} {amber :ALA}"
- some useful masks:
    - no solvent "not {solvent}"
    - just the protein "protein"
    - heavy protein atoms "diff {protein} {amber @/H}"
    - no hydrogens "not {amber @/H}"
    - pocket and ligand "pocket ligand"
    - specified ligand with 10 closest neighbours "resname UNL near UNL 10"
