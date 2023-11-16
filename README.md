
# EasyTrajH5

Trajectory management for mdtraj H5 files

## Installation

    pip install easytrajh5

## Quick Guide

We use the MD file format defined in the mdtraj library. We
have created `EasyTrajH5` which is a drop-in replacement
for `mdtraj.H5TrajectryFile` object with extra functionality.

```python
from easytrajh5.traj import EasyTrajH5File
 
h5 = EasyTrajH5File('traj.h5')
traj = h5.read_as_traj()
```

Load individual frames

```python
last_frame_traj = h5.read_frame_as_traj(-1)
```

As we use the `h5py` library, we can use efficient
fancy indexing to load just certain atoms:

```python
atom_indices = [100, 115, 116]
three_atom_traj = h5.read_as_traj(atom_indices=atom_indices)
```

We provide an atom mask selection function using a 
powerful new atom selection language (described in detail below):

```python
from easytrajh5.traj import EasyTrajH5File
 
mask = "intersect {mdtraj name CA} {protein}"
ca_trace_traj = EasyTrajH5File('traj.h5', atom_mask=mask).read_as_traj()
```

## Use as H5

There are convenience functions to insert different types
of data. 

To save/load strings:

```python
h5.set_str_dataset('my_string', 'a string')
new_str = h5.get_str_dataset('my_string')
```

To save/load json:
```python
h5.set_json_dataset('my_obj', {"a", "b"})
new_obj = h5.get_json_dataset('my_obj')
```
To insert/extract binary files:

```python
h5.insert_file_to_dataset('blob', 'blob.bin')
h5.extract_file_from_dataset('blob', 'new_blob.bin')
```

## EasyTrajH5 Atom Selection Language

Why another atom selection language since AMBER and MDTRAJ provides
one?. There are two main reasons. First, we wanted user-defined 
residue selections. These are stored in `easytrajh5/data/select.yaml`. 
Edit this file to create any new residue selections.

Second, we wanted to fix residue selection. The problem
is that AMBER uses residue numbering (':3,5,10-12') defined in the PDB file 
and not 0-based residue indexing. This means that in PDB files with multiple
chains, the residue number is not unique. MDTRAJ on the other hand, uses 
0-based indexing, but only allows you to use ranges ('resi 10 to 15'). We've
combined these ideas to provide our new flexible 0-based residue indexing 
'resi 3,5,10-12,100-150,300'.

We also allow you to easily drop in to AMBER and MDTRAJ simply by 
using the 'amber' and 'mdtraj' keywords. When combined with set 
operations, everything is now at your disposal.

Some useful masks:

- no solvent: "not {solvent}"
- just the protein: "protein"
- ligand and specific residues: "ligand resi 5,1,22-200"
- heavy protein atoms: "diff {protein} {amber @/H}"
- no hydrogens: "not {amber @/H}"
- ligand and 6 closest residues: "pocket ligand"
- specified ligand with 10 closest neighbours: "resname UNL near UNL 10"

#### Our user-defined keywords and special operator keywords

- accepts (in any order): 'ligand', 'protein', 'water', 'lipid', 'salt',
  'solvent', 'lipid', 'nucleic'
- If more than one keyword is specified, it is assumed they are joined with "or"
  operation (i.e. 'ligand protein' will return both ligand and protein atom indices).
- 'ligand' will find the residues 'LIG', 'UNL', 'UNK', or
   whatever is in 'ligand' in the h5 'easytrajh5/data/select.yaml'
- more can be defined in `easytrajh5/data/select.yaml`

Special operator keywords:

- 'pocket' will find the closest 6 residues to the 'ligand' group.
- 'near' will require a following resname, with an optional integer, e.g.:
    'near ATP'
    'near ATP 5'
- 'resname' identifies a single residue type
    'resname LEU'
- 'resi' for 0-indexed residue selections
    "resi 0,10-13" - selects atoms in the first and 11th to 14th residues
- 'atom' for 0-indexed atoms selections
    "atom 0,55,43,101-105" - selects the first, 56th, 44th, 102 to 106th atom

#### AMBER-style atom selection

- https://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax
- "amber :ALA,LYS" - selects all alanine and lysine residues

#### MDTraj-style atom selection 
- https://mdtraj.org/1.9.4/atom_selection.html
- "mdtraj protein and water" - selects protein and water

#### Set operations

Selections can be combined with set operators ("not", "intersect", "merge", "diff"):

- "intersect {not {amber :ALA}} {protein}"
- "diff {protein} {not {amber :ALA}}"
- "not {resname LEU}"
- "merge {near BSM 8} {amber :ALA}"


## Command-line utility `easyh5`

We have a command line `easyh5` to interrogate any `h5` file. To get a
schema of the dataset layout and attributes:

```bash
> easyh5 schema traj.h5
# {
# │   'datasets': [
# ....
# │   │   {
# │   │   │   'key': 'coordinates',
# │   │   │   'shape': [200, 3340, 3],
# │   │   │   'chunks': [3, 3340, 3],
# │   │   │   'is_extensible': True,
# │   │   │   'frame_shape': [3340, 3],
# │   │   │   'n_frame': 200,
# │   │   │   'dtype': 'float32',
# │   │   │   'attr': {'CLASS': 'EARRAY', 'EXTDIM': 0, 'TITLE': None, 'VERSION': '1.1', 'units': 'nanometers'}
# │   │   },
# ...
# │   │   {
# │   │   │   'key': 'topology',
# │   │   │   'shape': [1],
# │   │   │   'dtype': 'string(217329)',
# │   │   │   'attr': {'CLASS': 'ARRAY', 'FLAVOR': 'python', 'TITLE': None, 'VERSION': '2.4'}
# │   │   }
# │   ],
# │   'attr': {
# │   │   'CLASS': 'GROUP',
# │   │   'FILTERS': 65793,
# │   │   'PYTABLES_FORMAT_VERSION': '2.1',
# │   │   'TITLE': None,
# │   │   'VERSION': '1.0',
# │   │   'application': 'MDTraj',
# │   │   'conventionVersion': '1.1',
# │   │   'conventions': 'Pande',
# │   │   'program': 'MDTraj',
# │   │   'programVersion': '1.9.7',
# │   │   'title': 'title'
# │   }
# }
```

Or as a quick summary table:

```bash
> easyh5 size examples/trajectory.h5 
# 
#     examples/trajectory.h5     
# ┏━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┓
# ┃ dataset         ┃ size (MB) ┃
# ┡━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━┩
# │ cell_angles     │     <1 KB │
# │ cell_lengths    │     <1 KB │
# │ coordinates     │   7.64 MB │
# │ dry_atoms       │   0.03 MB │
# │ dry_topology    │   0.18 MB │
# │ kineticEnergy   │     <1 KB │
# │ potentialEnergy │     <1 KB │
# │ temperature     │     <1 KB │
# │ time            │     <1 KB │
# │ topology        │   0.21 MB │
# │ total           │   8.06 MB │
# └─────────────────┴───────────┘
```


## Miscellaneous utility functions

In `easytrajh5.struct`, there are some common transforms between MD prep 
and MD analysis packages:

```python
def dump_parmed(pmd: parmed.Structure, fname: str): 
def load_parmed(fname: str) -> parmed.Structure:
def get_parmed_from_pdb(pdb: str) -> parmed.Structure:
def get_parmed_from_parmed_or_pdb(pdb_or_parmed: str) -> parmed.Structure:
def get_parmed_from_mdtraj(traj: mdtraj.Trajectory, i_frame=0) -> parmed.Structure:
def get_parmed_from_openmm(openmm_topology, openmm_positions=None) -> parmed.Structure:
def get_mdtraj_from_parmed(pmd: parmed.Structure) -> mdtraj.Trajectory:
def get_mdtraj_from_openmm(openmm_topology, openmm_positions) -> mdtraj.Trajectory:
```

In `easytrajh5.quantity` we have some useful transforms to handle those
pesky unit objects from openmm. These transforms are used in our yaml and
json convenience functions

```python
from easytrajh5 import quantity
from parmed import unit

x = 5 * unit.nanosecond
d = quantity.get_dict_from_quantity(x)
# {
#│   'type': 'quantity',
#│   'value': 5,
#│   'unit': 'nanosecond',
#│   'unit_repr': 'Unit({BaseUnit(base_dim=BaseDimension("time"), name="nanosecond", symbol="ns"): 1.0})'
#}
y = quantity.get_quantity_from_dict(d)
# Quantity(value=5, unit=nanosecond)
```

