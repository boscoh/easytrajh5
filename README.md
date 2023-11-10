
# EasyTrajH5

Library for trajectory management with mdtraj H5 files

- EasyH5
    - simplified H5 class 
    - string methods
    - json methods
    - file-imports methods
    - easy dataset management
    - easy access to attr
    - schema dictionary
- EasyTrajH5File
    - drop-in replacement for mdtraj.H5TrajecotryFile
    - subclasses EasyH5, allowing access to underlying h5 methods
    - allows easy overrid to the h5py remote client
- select_mask
    - powerful atom selection language
    - proper set buildup using 'not', 'diff', 'merge', 'intersect'
    - allows amber and mdtraj style selections
    - proper 0-based indexing using flexible numbering list and ranges
    - user defined residue classes
- conveient transforms between mdtraj, parmed and openmm 

