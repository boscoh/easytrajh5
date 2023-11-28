#!/usr/bin/env python
import logging
import os

import chevron
import mdtraj
from path import Path

from easytrajh5.binary import get_binary
from easytrajh5.select import select_mask, slice_parmed_with_mask
from easytrajh5.struct import get_parmed_from_parmed_or_pdb
from easytrajh5.struct import load_parmed, slice_parmed, get_mdtraj_topology_from_pmd
from easytrajh5.traj import convert_h5_to_dcd_and_pdb

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

color1 = "green"
color2 = "purple"

colorid1 = 7
colorid2 = 11


def get_parmed_and_pdb(filename, is_solvent=True):
    """
    :param filename: str - either a .parmed or .pdb file
    :return: (parmed.Structure, str) - second str is the filename of a PDB file
    """
    filename = Path(filename)
    ext = filename.ext.lower()

    if ext not in [".pdb", ".parmed"]:
        raise ValueError(
            f"Can't process {filename} of type {ext}, only .pdb and .parmed"
        )

    pmd = get_parmed_from_parmed_or_pdb(filename)

    if ext == ".pdb" and is_solvent:
        return pmd, filename

    if not is_solvent:
        pmd = slice_parmed_with_mask(pmd, "not {solvent}")

    pdb = filename.with_suffix(".parmed.pdb")
    pmd.save(pdb, overwrite=True)

    return pmd, pdb


pymol_template = """
set defer_builds_mode, 3
{{load_script}}
hide
select backbone, name c+o+n+h+oxt and not solvent
deselect

{{#is_stick}}
show sticks
{{/is_stick}}

{{^is_stick}}
show sticks
hide sticks, backbone
hide sticks, hydro
show sticks, name CA
set cartoon_flat_sheets, 0
set cartoon_loop_radius, 0.4
show cartoon
{{/is_stick}}

{{#is_sphere}}
show spheres
{{/is_sphere}}

util.cba(29,"{{name}}")

{{#ranks}}
select rank {{ranks}}
color {{color1}}, sele
deselect
{{/ranks}}

{{#ranks2}}
select rank {{ranks2}}
color {{color2}}, sele
deselect
{{/ranks2}}
"""


def pymol(filename, mask, mask2, sphere, stick, keep_solvent):
    name = Path(filename).name
    if filename.endswith(".h5"):
        dcd, pdb = convert_h5_to_dcd_and_pdb(filename, keep_solvent)
        structure, pdb = get_parmed_and_pdb(pdb)
        load_script = f"load {pdb}, {name}\n"
        load_script += f"load_traj {dcd}, {name}\n"
    else:
        structure, pdb = get_parmed_and_pdb(filename, keep_solvent)
        load_script = f"load {pdb}, {name}\n"

    data = {
        "load_script": load_script,
        "is_stick": stick,
        "is_sphere": sphere,
        "name": name,
        "color1": color1,
        "color2": color2,
    }
    if mask:
        i_atoms = select_mask(structure, mask)
        data["ranks"] = "+".join(map(str, i_atoms))
    if mask2:
        i_atoms = select_mask(structure, mask2)
        data["ranks2"] = "+".join(map(str, i_atoms))

    logger.info("Starting pymol")
    pymol_script_fname = Path(filename).with_suffix(".pml")
    with open(pymol_script_fname, "w") as f:
        f.write(chevron.render(template=pymol_template, data=data))
    os.system(f"{get_binary('pymol')} {pymol_script_fname}")


chimera_template = """
~display

{{load_script}}

repr stick
~display H
setattr m stickScale 2

{{#is_sphere}}
~ribbon
repr cpk
display
{{/is_sphere}}

{{#is_stick}}
repr stick
~ribbon
display
{{/is_stick}}

{{#select_str}}
color {{color1}},a,r {{select_str}}
{{/select_str}}

{{#select_str2}}
color {{color2}},a,r {{select_str2}}
{{/select_str2}}

start Side View

"""


def chimera(filename, parmed, mask, mask2, sphere, stick, keep_solvent):
    if filename.endswith(".h5"):
        dcd, pdb = convert_h5_to_dcd_and_pdb(filename, keep_solvent)

        if not parmed:
            raise ValueError("Error: -p/--parmed must be set")
        parmed = Path(parmed)
        amber_parm7 = parmed.with_suffix(".parm7")
        pmd = load_parmed(parmed)
        if not keep_solvent:
            pmd = slice_parmed(pmd, select_mask(pmd, "not {solvent}"))
        pmd.save(amber_parm7, overwrite=True)

        chimera_load = Path(filename).with_suffix(".chimera.load")
        with open(chimera_load, "w") as f:
            f.write("\n".join(["namdprmtopdcd", amber_parm7, dcd, ""]))
        load_script = f"open movie:{chimera_load}"
    else:
        pmd, pdb = get_parmed_and_pdb(filename, keep_solvent)
        load_script = f"open {pdb}"

    data = {
        "load_script": load_script,
        "is_sphere": sphere,
        "is_stick": stick,
        "color1": color1,
        "color2": color2,
    }
    if mask:
        i_atoms = select_mask(pmd, mask)
        atoms = []
        top = get_mdtraj_topology_from_pmd(pmd)
        for r in top.residues:
            names = [a.name for a in r.atoms if a.index in i_atoms]
            if len(names):
                atoms.append(f":{r.index + 1}@{','.join(names)}")
        data["select_str"] = " ".join(atoms)
    if mask2:
        i_atoms = select_mask(pmd, mask2)
        atoms = []
        top = get_mdtraj_topology_from_pmd(pmd)
        for r in top.residues:
            names = [a.name for a in r.atoms if a.index in i_atoms]
            if len(names):
                atoms.append(f":{r.index + 1}@{','.join(names)}")
        data["select_str2"] = " ".join(atoms)

    logger.info("Starting chimera")
    chimera_cmd_fname = Path(filename).with_suffix(".chimera.cmd")
    with open(chimera_cmd_fname, "w") as f:
        f.write(chevron.render(template=chimera_template, data=data))
    os.system(f"{get_binary('chimera')} {chimera_cmd_fname}")


vmd_template = """
# turn on lights 0 and 1
light 0 on
light 1 on
light 2 off
light 3 off

# position the stage and axes
axes location off
stage location off

# position and turn on menus
menu main on
menu graphics on
menu main move 700 60
menu graphics move 700 340

display projection orthographic
display depthcue on
display cuestart 1.5
display cueend 2.75
display cuemode Linear

mol new "{{top}}" type {{top_type}}
mol delrep 0 0
{{#trj}}
mol addfile "{{trj}}" type {{trj_type}}
{{/trj}}


{{#is_sphere}}
mol selection {all}
mol representation VDW
mol addrep top

{{#ranks}}
mol selection {{ranks}}
mol representation VDW
mol color colorid {{colorid1}}
mol addrep top
{{/ranks}}

{{#ranks2}}
mol selection {{ranks2}}
mol representation VDW
mol color colorid {{colorid2}}
mol addrep top
{{/ranks2}}

{{/is_sphere}}


{{^is_sphere}}

mol selection {backbone}
mol representation NewCartoon 0.40 20.00 2.50 0
mol addrep top
mol selection {(not (protein and backbone or water) or name CA) and not name "[0-9]?H.*"}
mol representation CPK 1.4 2 8 6
mol addrep top

{{#ranks}}
mol selection {{ranks}}
mol representation CPK 1.4 2 8 6
mol color colorid {{colorid1}}
mol addrep top
{{/ranks}}

{{#ranks2}}
mol selection {{ranks2}}
mol representation CPK 1.4 2 8 6
mol color colorid {{colorid2}}
mol addrep top
{{/ranks2}}

{{/is_sphere}}
"""


def vmd(filename, mask, mask2, sphere, keep_solvent):
    if filename.endswith(".h5"):
        dcd, pdb = convert_h5_to_dcd_and_pdb(filename, keep_solvent)
        structure, pdb = get_parmed_and_pdb(pdb, keep_solvent)
    else:
        dcd = ""
        structure, pdb = get_parmed_and_pdb(filename, keep_solvent)

    data = {
        "top_type": "pdb",
        "trj_type": "dcd",
        "top": pdb,
        "trj": dcd,
        "is_sphere": sphere,
        "colorid1": colorid1,
        "colorid2": colorid2,
    }
    if mask:
        i_atoms = select_mask(structure, mask)
        data["ranks"] = "{index " + " ".join(map(str, i_atoms)) + "}"
    if mask2:
        i_atoms = select_mask(structure, mask2)
        data["ranks2"] = "{index " + " ".join(map(str, i_atoms)) + "}"

    logger.info("Starting vmd")
    vmd_script_fname = Path(filename).with_suffix(".vmd")
    with open(vmd_script_fname, "w") as f:
        f.write(chevron.render(template=vmd_template, data=data))
    os.system(f"{get_binary('vmd')} -e {vmd_script_fname}")
