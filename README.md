
# EasyTrajH5

Trajectory management with mdtraj H5 files

- EasyH5File
  - methods-based interface for H5
  - string, json, file-import methods
  - easy dataset/attr management
  - schema dictionary
- EasyTrajH5File
  - drop-in replacement for mdtraj.H5TrajecotryFile
  - subclasses EasyH5File
  - allows simple override to h5pyd
- select_mask
  - powerful atom selection language
  - proper set buildup using 'not', 'diff', 'merge', 'intersect'
  - allows amber and mdtraj style selections
  - proper 0-based indexing using flexible numbering list and ranges
  - user defined residue classes
- conveient transforms
  - mdtraj
  - parmed 
  - openmm 
- cli for interrogating h5 files


### EasyH5

Accessing H5 file for basic data needs