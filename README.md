
# EasyTrajH5

Trajectory management for mdtraj H5 files with atom selection language
and efficient data operations via the h5py library.

## Installation

    pip install easytrajh5

## Quick Guide

Our main file object `EasyTrajH5` is a drop-in replacement
for `mdtraj.H5TrajectryFile`:

```python
from easytrajh5.traj import EasyTrajH5File
 
h5 = EasyTrajH5File('traj.h5')
traj = h5.read_as_traj()
```
This loads the data progressively in chunks, allowing online streaming 
in advanced usage. 

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

We provide atom selection using a new selection language (described in detail below).
This is particular efficient as it only loads the atoms you want, without 
requiring the entire trajectory to be loaded into memory:

```python
from easytrajh5.traj import EasyTrajH5File
 
mask = "intersect {mdtraj name CA} {protein}"
ca_trace_traj = EasyTrajH5File('traj.h5', atom_mask=mask).read_as_traj()
```

Drop in replacement for `mdtraj.reporters.HDF5Reporter` in openmm
that uses `EasyTrajH5File`: 

```python
from easytrajh5.traj import EasyTrajH5Reporter
```


## Atom Selection Language

Why another atom selection language (we have AMBER and MDTRAJ)?
Two main reasons. 

First, we wanted user-defined 
residue selections. These are stored in `easytrajh5/data/select.yaml`. 
Edit this file to create any new residue selections.

Second, we wanted to fix residue selection. The problem
is that AMBER uses residue numbering (`:3,5,10-12`) defined in the PDB file 
and not 0-based residue indexing. This means that in PDB files with multiple
chains, the residue number is not unique. MDTRAJ on the other hand, uses 
0-based indexing, but only allows you to use ranges (`resi 10 to 15`). 

We've combined these ideas to provide our new flexible 0-based residue indexing 
`resi 3,5,10-12,100-150,300`.

We also allow you to easily drop in to AMBER and MDTRAJ simply by 
using the `amber` and `mdtraj` keywords. When combined with set 
operations, everything is now at your disposal.

Some useful masks:

- no solvent: `not {solvent}`
- just the protein: `protein`
- ligand and specific residues: `ligand resi 5,1,22-200`
- heavy protein atoms: `diff {protein} {amber @/H}`
- no hydrogens: `not {amber @/H}`
- ligand and 6 closest residues: `pocket ligand`
- specified ligand with 10 closest neighbours: `resname UNL near UNL 10`

#### User-defined and operator keywords

If more than one keyword is specified, it is assumed they are joined with "or"
operation (i.e. `ligand protein` will return both ligand and protein atom indices).

This default keywords are:
- `ligand`, `protein`, `water`, `lipid`, `salt`, `solvent`, `lipid`, `nucleic`
- as defined in `easytrajh5/data/select.yaml`
- `ligand` will find the residues `LIG`, `UNL`, `UNK`

Special operator keywords:

- `pocket` will find the closest 6 residues to the `ligand` group.
- `near` will require a following resname, with an optional integer, e.g.:
    `near ATP`
    `near ATP 5`
- `resname` identifies a single residue type
    `resname LEU`
- `resi` for 0-indexed residue selections
    `resi 0,10-13` - selects atoms in the first and 11th to 14th residues
- `atom` for 0-indexed atoms selections
    `atom 0,55,43,101-105` - selects the first, 56th, 44th, 102 to 106th atom

#### AMBER-style atom selection

- https://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax
- `amber :ALA,LYS` - selects all alanine and lysine residues

#### MDTraj-style atom selection 
- https://mdtraj.org/1.9.4/atom_selection.html
- `mdtraj protein and water` - selects protein and water

#### Set operations

Selections can be combined with set operators: `not`, `intersect`, `merge`, `diff`:

- `intersect {not {amber :ALA}} {protein}`
- `diff {protein} {not {amber :ALA}}`
- `not {resname LEU}`
- `merge {near BSM 8} {amber :ALA}`

#### Use in python

In your python code, there is a `select_mask` fn that operates on `parmed.Structure`
objects:

```python
from easytrajh5.traj import EasyTrajH5File
from easytrajh5.select import select_mask
from easytrajh5.struct import slice_parmed

pmd = EasyTrajH5File("traj.h5").get_topology_parmed()
i_atoms = select_mask(pmd, "not {solvent}")
sliced_pmd = slice_parmed(pmd, i_atoms)
```

Some common conversions and loaders in `easytrajh5.struct` for `parmed.Structure` and
`mdtraj.Trajectory` objects:

```python
import parmed, mdtraj

def dump_parmed(pmd: parmed.Structure, fname: str): 
def load_parmed(fname: str) -> parmed.Structure:
def get_parmed_from_pdb(pdb: str) -> parmed.Structure:
def get_parmed_from_parmed_or_pdb(pdb_or_parmed: str) -> parmed.Structure:
def get_parmed_from_mdtraj(traj: mdtraj.Trajectory, i_frame=0) -> parmed.Structure:
def get_parmed_from_openmm(openmm_topology, openmm_positions=None) -> parmed.Structure:
def get_mdtraj_from_parmed(pmd: parmed.Structure) -> mdtraj.Trajectory:
def get_mdtraj_from_openmm(openmm_topology, openmm_positions) -> mdtraj.Trajectory:
```

## Use as H5

There are convenience functions to insert different types
of data. 

To save/load strings:

```python
h5.set_str_dataset('my_string', 'a string')
h5.flush()
new_str = h5.get_str_dataset('my_string')
```

To save/load json:
```python
h5.set_json_dataset('my_obj', {"a", "b"})
h5.flush()
new_obj = h5.get_json_dataset('my_obj')
```
To insert/extract binary files:

```python
h5.insert_file_to_dataset('blob', 'blob.bin')
h5.flush()
h5.extract_file_from_dataset('blob', 'new_blob.bin')
```

We can get information about the h5 file:

```python
schema_json = h5.get_schema()
dataset_keys = h5.get_dataset_keys()
attr_keys = h5.get_attr_keys()
```

We can extract data

```python
dataset = h5.get_dataset("coordinates")
value_list = dataset[:]
last_value = dataset[-1]

# if the attrs are set
value = h5.get_attr('user')
```

Convenience function to append values to an `h5` file without
worrying about file or dataset creation:

```python
from easytrajh5.h5 import dump_value_to_h5, EasyH5File

dump_value_to_h5('new.h5', [1,2], 'my_data_set')
dump_value_to_h5('new.h5', [3,4], 'my_data_set')
dump_value_to_h5('new.h5', [5,7], 'my_data_set')

return_values = EasyH5File('new.h5').get_dataset("my_data_set")[:]
# [[1,2], [3,4], [5,6]]
```

## Command-line utility `easyh5`

`easyh5` provides a bunch of useful cli subcommands to interrogate `h5` and related files:

```
Usage: easyh5 [OPTIONS] COMMAND [ARGS]...

  h5: preprocessing and analysis tools

Options:
  --help  Show this message and exit.

Commands:
  dataset        Examine contents of h5
  insert-parmed  Insert parmed into dataset:parmed of an H5
  mask           Explore residues/atoms of H5/PDB/PARMED using mask
  merge          Merge a list of H5 files
  parmed         Extract parmed from dataset:parmed of an H5 with...
  pdb            Extract PDB of a frame of an H5
  schema         Examine layout of H5
  show-chimera   Use CHIMERA to show H5/PDB/PARMED with mask, needs PARMED
  show-pymol     Use PYMOL to show H5/PDB/PARMED with mask
  show-vmd       Use VMD to show H5/PDB/PARMED with mask
```

To get a schema of the dataset layout and attributes:

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
#
# ...
#
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
> easyh5 dataset examples/trajectory.h5 
# Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
# 
#                   sims/high_bf/trajectory.h5                  
#                                                               
#   dataset           shape              dtype       size (MB)  
#  ──────────────────────────────────────────────────────────── 
#   cell_angles       (1500, 3)          float32       0.02 MB  
#   cell_lengths      (1500, 3)          float32       0.02 MB  
#   coordinates       (1500, 25767, 3)   float32     442.32 MB  
#   kineticEnergy     (1500,)            float32         <1 KB  
#   potentialEnergy   (1500,)            float32         <1 KB  
#   temperature       (1500,)            float32         <1 KB  
#   time              (1500,)            float32         <1 KB  
#   topology          (1,)               |S2083249     1.99 MB  
#                                                               
#   total                                            444.36 MB  
```

To get an overview of a dataset:

```bash
> easyh5 dataset examples/trajectory.h5 coordinates
# Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
# 
#   examples/trajectory
#      dataset=coordinates
#      shape=(1500, 25767, 3)
# 
# [[[1.291678   7.558739   1.5199517 ]
#   [1.368739   7.5888386  1.4620152 ]
#   [1.2175218  7.6268845  1.5275735 ]
#   ...
#   [2.375777   0.09478953 4.0356894 ]
#   [3.107005   3.3255231  2.8464174 ]
#   [3.0329072  3.9307644  1.3600407 ]]
# 
#  ...
# 
#  [[2.9693408  7.1466036  1.4656581 ]
#   [2.9327238  7.198606   1.3871984 ]
#   [3.0665123  7.171176   1.4781022 ]
#   ...
#   [4.944392   0.56028575 4.301907  ]
#   [2.6180382  0.3969128  1.4842175 ]
#   [3.281546   4.9666233  2.4855924 ]]]
```

Or to focus on a selected frames, use a numbered lis:

```bash
> easyh5 dataset examples/trajectory.h5 coordinates 1,3,4-10
# Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
# 
#   sims/high_bf/trajectory.h5
#      dataset=coordinates
#      shape=(1500, 25767, 3)
# 
# frames(1,3,4-10)=
# [[[1.2958181  7.5481067  1.5513833 ]
#   [1.2766361  7.4586782  1.5085387 ]
#   [1.3654946  7.58469    1.4880756 ]
#   ...
#   [2.0149727  0.20826703 3.712016  ]
#   [3.3603299  3.6615734  2.6487541 ]
#   [3.1595583  4.0199933  1.509442  ]]
# 
#  ...
# 
#  [[1.2550778  7.4836254  1.5989571 ]
#   [1.278228   7.403505   1.5419852 ]
#   [1.2919694  7.571082   1.5644412 ]
#   ...
#   [2.6036768  5.9193387  3.7148886 ]
#   [4.2752028  3.8813443  2.6205144 ]
#   [2.343824   3.9689744  0.05281828]]]
# 
```

To check atom selections of the protein:

```bash
> easyh5 mask sims/high_bf/trajectory.h5 "amber :PRO" --res
# Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
# EasyTrajH5File: fname='sims/high_bf/trajectory.h5' mode='a' atom_mask='' is_dry_cache=False
# open connection:started...
# open connection:finished in <1ms
# loading topology:started...
# loading topology:finished in 134ms
# select_mask "amber :PRO" -> 112 atoms, 8 residues
# <Residue PRO[7]; chain=1>
# <Residue PRO[12]; chain=1>
# <Residue PRO[42]; chain=1>
# <Residue PRO[48]; chain=1>
# <Residue PRO[50]; chain=1>
# <Residue PRO[87]; chain=1>
# <Residue PRO[129]; chain=1>
# <Residue PRO[167]; chain=1>
```

To extract that as PDB:

```bash
> easyh5 mask sims/high_bf/trajectory.h5 "amber :PRO" --pdb pro.pdb
```

There are three sub-commands that help visualize selections in standard viewers:

- `easyh5 show-pymol <PDB> <MASK1> <MASK2>`
- `easyh5 show-vmd <PDB|PARMED|H5> <MASK1> <MASK2>`
- `easyh5 show-chimera <PDB|PARMED|H5> <MASK1> <MASK2>`

It will open the structure or trajectory in the corresponding viewers with the first selection
colored in green, and the second selection in pink.

A configuration file in your systems config directory `rseed.binary.yaml` will be created
that list the full path name of PYMOL/VMD/CHIMERA. Change this if your copy of the viewer 
is in a different location.


## Miscellaneous utility 

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

