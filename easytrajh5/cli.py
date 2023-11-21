#!/usr/bin/env python
import logging

import click
import mdtraj
from path import Path

from easytrajh5.h5 import EasyH5File
from easytrajh5.manager import TrajectoryManager
from easytrajh5.select import get_n_residue_of_mask, select_mask
from easytrajh5.select import parse_number_list
from easytrajh5.struct import slice_parmed, get_parmed_from_parmed_or_pdb
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
    """Examine layout of H5"""
    EasyH5File(h5).print_schema()


@h5.command()
@click.argument("h5")
@click.argument("dataset", default=None, required=False)
@click.argument("frames", default=None, required=False)
def dataset(h5, dataset, frames):
    """Examine contents of h5"""
    f = EasyH5File(h5)
    if dataset is None:
        f.print_dataset_table(h5)
    else:
        print(f"\n  {h5}")
        f.print_dataset(dataset, frames)


@h5.command()
@click.argument("h5")
@click.argument("dataset", required=False)
def json(h5, dataset):
    """
    Get JSON configs of H5
    """
    EasyH5File(h5).print_json(dataset)


@h5.command()
@click.argument("h5", default="trajectory.h5")
@click.option("--mask", default=None, help="atom selection", show_default=True)
@click.option("--i", default=0, help="frame", show_default=True)
def pdb(h5_trajectory, mask, i):
    """
    Extract PDB of a frame of an H5
    """
    pmd = EasyTrajH5File(h5_trajectory, atom_mask=mask).get_parmed(i_frame=i)
    pdb = Path(h5_trajectory).with_suffix(".pdb")
    pmd.save(pdb, overwrite=True)
    print(f"Wrote {pdb=} {i=} {mask=}")


@h5.command()
@click.argument("h5-pdb-parmed")
@click.argument("mask", default=None, required=False)
@click.option("--pdb", help="Save to PDB")
@click.option("--atom", flag_value=True, help="List all atoms")
@click.option("--res", flag_value=True, help="List all residues")
@click.option("--i", default=0, help="frame", show_default=True)
def mask(h5_pdb_parmed, mask, pdb, atom, res, i):
    """
    Explore residues/atoms of H5/PDB/PARMED using mask
    """
    filename = Path(h5_pdb_parmed)
    ext = filename.ext.lower()

    if ext not in [".h5", ".pdb", ".parmed"]:
        print(f"Can't recognize file extension {h5_pdb_parmed}")
        return

    if ext == ".h5":
        pmd = EasyTrajH5File(h5_pdb_parmed).get_parmed(i_frame=i)
    else:
        pmd = get_parmed_from_parmed_or_pdb(h5_pdb_parmed)

    if mask is None:
        get_n_residue_of_mask(pmd, "protein")
        get_n_residue_of_mask(pmd, "ligand")
        get_n_residue_of_mask(pmd, "solvent")
        get_n_residue_of_mask(pmd, "not {merge {protein} {solvent} {ligand}}")
        return

    i_atoms = select_mask(pmd, mask, is_fail_on_empty=False)
    if not len(i_atoms):
        print("Could select any atoms")
        return

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

    if not pdb:
        return

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
    """Merge a list of H5 files"""
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
