#!/usr/bin/env python
import logging
from pathlib import Path

import click
import mdtraj
import mdtraj as md
import numpy

from easytrajh5.fs import ensure_dir
from easytrajh5.h5 import EasyH5File, print_schema, print_size, print_json
from easytrajh5.manager import TrajectoryManager
from easytrajh5.select import get_n_residue_of_mask
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
    fname = str(Path(prefix).with_suffix(".h5"))
    print(f"Merged {h5_list} --> {fname}")
    frames.save_hdf5(fname)


@h5.command()
@click.argument("h5")
def schema(h5):
    """Examine contents of h5"""
    print_schema(EasyH5File(h5))


@h5.command()
@click.argument("h5")
def size(h5):
    """Examine contents of h5"""
    print_size(EasyH5File(h5), h5)


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
@click.option("--mask", "-m", default="protein")
@click.option("--output-dir", "-o", default=".")
def contacts(h5_trajectory, mask, output_dir):
    """
    Evaluates the residues selected with MASK in H5_TRAJECTORY to determine
    the average and standard deviation inter residue contact distance.
    Each is saved in the specified OUTPUT-DIR to their own NxN numpy array where N
    is the number of residues selected.

    \b
    NOTE: Distance calculated is between closest non H atoms between residues

    \b
    outputs:
      contact_avg.npy:  NxN numpy array of average contact
                        distance between residues
      contact_std.npy:  NxN numpy array of standard deviation contact
                        distance between residues
    """
    frames = EasyTrajH5File(h5_trajectory, atom_mask=mask).read_as_traj()

    anchors = [mol for mol in frames.topology.find_molecules()]
    anchors_united = set.union(*anchors)
    frames.image_molecules(inplace=True, anchor_molecules=[anchors_united])
    cont = md.compute_contacts(frames, contacts="all", ignore_nonprotein=False)

    sqr = md.geometry.squareform(cont[0], cont[1])

    ensure_dir(output_dir)
    fname_avg = Path(output_dir) / "contact_avg.npy"
    fname_std = Path(output_dir) / "contact_std.npy"
    with open(fname_avg, "wb") as f:
        numpy.save(f, numpy.mean(sqr, axis=0))
    with open(fname_std, "wb") as f:
        numpy.save(f, numpy.std(sqr, axis=0))


@h5.command()
@click.argument("h5-trajectory", default="trajectory.h5")
def residues(h5_trajectory):
    """
    Identify the types of residues in the mdtraj h5 file
    """
    pmd = EasyTrajH5File(h5_trajectory).get_parmed()
    print(f"residues in {h5_trajectory}")
    print(f'protein = {get_n_residue_of_mask(pmd, "protein")}')
    print(f'ligand = {get_n_residue_of_mask(pmd, "ligand")}')
    print(f'solvent = {get_n_residue_of_mask(pmd, "solvent")}')
    print(
        f'other = {get_n_residue_of_mask(pmd, "not {merge {protein} {solvent} {ligand}}")}'
    )


if __name__ == "__main__":
    h5()
