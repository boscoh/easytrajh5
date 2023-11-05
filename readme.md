# EasyTrajH5

Library to open-source for trajectory management using mdtraj h5

- EasyH5
    - simplified H5 class EasyH5
    - string methods
    - json methods
    - file-imports methods
    - simple dataset creation/append
- EasyTrajH5File
    - drop-in replacement for mdtraj.H5TrajecotryFile
    - subclasses EasyH5, allowing access to h5 methods
    - uses h5py, allows easy override for h5pyd
- select_mask
    - powerful atom selection language
- transforms between mdtraj, parmed and openmm topology/position objects
- pdb editing functions

Dependencies:

- addict==2.4.0
- deepdiff==6.3.1
- h5py==3.8.0
- mdtraj==1.9.7
- numpy==1.23.3
- orjson==3.9.10
- ParmEd==3.4.3
- pydash==7.0.6
- rich==13.6.0
- rseed.egg==info
- ruamel.base==1.0.0