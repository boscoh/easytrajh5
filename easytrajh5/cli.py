#!/usr/bin/env python
import logging

import click
import mdtraj
from path import Path

from easytrajh5.h5 import EasyH5File, print_schema, print_size, print_json
from easytrajh5.manager import TrajectoryManager
from easytrajh5.select import get_n_residue_of_mask, select_mask
from easytrajh5.struct import slice_parmed
from easytrajh5.traj import EasyTrajH5File

logger = logging.getLogger(__name__)

logging.basicConfig(format="%(message)s", level=logging.INFO)


@click.group()
def h5():
    """
    h5: preprocessing and analysis tools
    """
    pass


@h5.command()
@click.argument("h5")
def schema(h5):
    """Examine layout of h5"""
    print_schema(EasyH5File(h5))


@h5.command()
@click.argument("h5")
def datasets(h5):
    """Examine contents of h5"""
    print_size(EasyH5File(h5), h5)


@h5.command()
@click.argument("h5")
@click.argument("dataset")
@click.argument("frames", default=None, required=False)
def peek(h5, dataset, frames):
    """Examine contents of h5"""
    from easytrajh5.select import parse_number_list

    f = EasyH5File(h5)
    d = f.get_dataset(dataset)
    print(f"Data at {h5}/{dataset}/shape{d.shape}")
    if frames is not None:
        i_frames = parse_number_list(frames)
        chunk = d[i_frames]
        print(f"Peek at {frames} shape{chunk.shape}")
    else:
        chunk = d[:]
    print(chunk)


@h5.command()
@click.argument("h5")
@click.argument("dataset", required=False)
def json(h5, dataset):
    """
    Get JSON configs associated with entry
    """
    print_json(EasyH5File(h5), dataset)


@h5.command()
@click.argument("h5-trajectory", default="trajectory.h5")
def content(h5_trajectory):
    """
    Identify the types of residues in the mdtraj h5 file
    """
    pmd = EasyTrajH5File(h5_trajectory).get_parmed()
    get_n_residue_of_mask(pmd, "protein")
    get_n_residue_of_mask(pmd, "ligand")
    get_n_residue_of_mask(pmd, "solvent")
    get_n_residue_of_mask(pmd, "not {merge {protein} {solvent} {ligand}}")


@h5.command()
@click.argument("h5-trajectory", default="trajectory.h5")
@click.option(
    "--mask",
    default="",
    help="atom selection",
    show_default=True,
)
@click.option(
    "--i",
    default=0,
    help="frame",
    show_default=True,
)
def pdb(h5_trajectory, mask, i):
    """
    Extract PDB of a frame
    """
    pmd = EasyTrajH5File(h5_trajectory, atom_mask=mask).get_parmed(i_frame=i)
    pdb = Path(h5_trajectory).with_suffix(".pdb")
    pmd.save(pdb, overwrite=True)
    print(f"Wrote {pdb=} {i=} {mask=}")


@h5.command()
@click.argument("h5-trajectory", default="trajectory.h5")
@click.argument(
    "mask",
    default="",
)
@click.option("--pdb", help="Pdb to save")
@click.option("--atom", flag_value=True, help="List all atoms")
@click.option("--res", flag_value=True, help="List all residues")
@click.option(
    "--i",
    default=0,
    help="frame",
    show_default=True,
)
def mask(h5_trajectory, mask, pdb, atom, res, i):
    """
    Explore atoms and residues using the atom selection language
    """
    pmd = EasyTrajH5File(h5_trajectory, atom_mask=mask).get_parmed(i_frame=i)
    i_atoms = select_mask(pmd, mask)
    pmd = slice_parmed(pmd, i_atoms)
    if res:
        residues = []
        for a in pmd.atoms:
            if a.residue not in residues:
                residues.append(a.residue)
        for residue in residues:
            print(residue)
    if atom:
        for a in pmd.atoms:
            print(a)
    if pdb:
        pdb = Path(pdb).with_suffix(".pdb")
        pmd.save(pdb, overwrite=True)
        print(f"Wrote {pdb}")


@h5.command()
@click.argument("h5_list", nargs=-1)
@click.option(
    "--prefix",
    default="merged",
    help="prefix for newly generated .h5",
    show_default=True,
)
@click.option(
    "--mask",
    default="amber @*",
    help="selection mask to specify atoms in newly generated .h5",
    show_default=True,
)
def merge(h5_list, prefix, mask):
    """Merge a list of h5 files into one PREFIX.h5, a subset of atoms that exist in all files
    can be specified with MASK"""
    traj_mananger = TrajectoryManager(paths=h5_list, atom_mask=f"{mask}")
    frames = []
    for t_id in traj_mananger.traj_file_by_i.keys():
        for f_id in range(0, traj_mananger.get_n_frame(t_id)):
            frames.append(traj_mananger.read_as_frame_traj((f_id, t_id)))
    frames = mdtraj.join(frames)
    fname = Path(prefix).with_suffix(".h5")
    print(f"Merged {h5_list} --> {fname}")
    frames.save_hdf5(fname)


if __name__ == "__main__":
    h5()
