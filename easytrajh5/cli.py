#!/usr/bin/env python
import logging

import click
import mdtraj
from path import Path

import easytrajh5.show
from easytrajh5.h5 import EasyH5File
from easytrajh5.manager import TrajectoryManager
from easytrajh5.select import get_n_residue_of_mask, select_mask
from easytrajh5.struct import slice_parmed, get_parmed_from_parmed_or_pdb, load_parmed
from easytrajh5.traj import EasyTrajH5File

logger = logging.getLogger(__name__)

logging.basicConfig(format="%(message)s", level=logging.INFO)


@click.group()
def h5():
    """
    h5: preprocessing and analysis tools
    """
    pass


@h5.command(no_args_is_help=True)
@click.argument("h5")
def schema(h5):
    """Examine layout of H5"""
    EasyH5File(h5).print_schema()


@h5.command(no_args_is_help=True)
@click.argument("h5")
@click.argument("dataset", default=None, required=False)
@click.option(
    "--frames", default=None, help="list of frames e.g. '1,2,3-10,11,200-211'"
)
@click.option("--json", is_flag=True, default=False, help="Format json")
def dataset(h5, dataset, frames, json):
    """Examine contents of h5"""
    f = EasyH5File(h5)
    title = h5
    if dataset is None:
        f.print_dataset_table(title)
    else:
        print(f"\n  {title}")
        f.print_dataset(dataset, frames=frames, is_json=json)


@h5.command(no_args_is_help=True)
@click.argument("h5", default="trajectory.h5")
@click.option("--mask", default=None, help="atom selection", show_default=True)
@click.option("--i", "-i", default=0, help="frame", show_default=True)
def pdb(h5_trajectory, mask, i):
    """
    Extract PDB of a frame of an H5
    """
    pdb = Path(h5_trajectory).with_suffix(".pdb")
    pmd = EasyTrajH5File(h5_trajectory, atom_mask=mask).get_topology_parmed(i_frame=i)
    pmd.save(pdb, overwrite=True)
    print(f"Wrote {pdb=} {i=} {mask=}")


@h5.command(no_args_is_help=True)
@click.argument("h5")
@click.argument("parmed")
@click.argument("dataset", default="parmed")
def insert_parmed(h5, parmed, dataset):
    """
    Insert parmed into dataset:parmed of an H5
    """
    if dataset is None:
        dataset = "parmed"
    in_pmd = load_parmed(parmed)
    n_atom_in_parmed = len(in_pmd.atoms)
    with EasyTrajH5File(h5) as traj_file:
        assert n_atom_in_parmed == traj_file.get_n_atom()
        traj_file.insert_file_to_dataset(dataset, parmed)
        logger.info(f"Inserted {parmed} to {h5}:parmed")


@h5.command(no_args_is_help=True)
@click.argument("h5")
@click.option("--mask", default=None, help="atom selection", show_default=True)
@click.option("--i", "-i", default=0, help="frame", show_default=True)
def parmed(h5, mask, i):
    """
    Extract parmed from dataset:parmed of an H5 with optional frame
    """
    pdb = Path(h5).with_suffix(".pdb")
    pmd = EasyTrajH5File(h5, atom_mask=mask).get_parmed_from_dataset(i_frame=i)
    pmd.save(pdb, overwrite=True)
    print(f"Wrote {pdb=} {i=} {mask=}")


@h5.command(no_args_is_help=True)
@click.argument("h5-pdb-parmed")
@click.argument("mask", default=None, required=False)
@click.option("--pdb", help="Save to PDB")
@click.option("--atom", flag_value=True, help="List all atoms")
@click.option("--res", flag_value=True, help="List all residues")
@click.option("--i", "-i", default=0, help="frame", show_default=True)
def mask(h5_pdb_parmed, mask, pdb, atom, res, i):
    """
    Explore residues/atoms of H5/PDB/PARMED using mask
    """
    filename = Path(h5_pdb_parmed)
    ext = filename.suffix.lower()

    if ext not in [".h5", ".pdb", ".parmed"]:
        print(f"Can't recognize file extension {h5_pdb_parmed}")
        return

    if ext == ".h5":
        pmd = EasyTrajH5File(h5_pdb_parmed).get_topology_parmed(i_frame=i)
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
        print("Couldn't select any atoms")
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


@h5.command(no_args_is_help=True)
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
    traj_mananger = TrajectoryManager(paths=h5_list, atom_mask=mask)
    frames = []
    for i_traj in traj_mananger.traj_file_by_i.keys():
        for i_frame in range(0, traj_mananger.get_n_frame(i_traj)):
            frames.append(traj_mananger.read_as_frame_traj((i_frame, i_traj)))
    full_traj = mdtraj.join(frames)
    full_traj_h5 = Path(prefix).with_suffix(".h5")
    print(f"Merged {h5_list} --> {full_traj_h5}")
    full_traj.save_hdf5(full_traj_h5)


@h5.command(no_args_is_help=True)
@click.argument("filename")
@click.argument("mask", required=False)
@click.argument("mask2", required=False)
@click.option("-s", "--sphere", is_flag=True, help="Show all atoms as spheres")
@click.option("-t", "--stick", is_flag=True, help="Show all bonds as stick")
@click.option(
    "-k",
    "--keep-solvent",
    is_flag=True,
    default=False,
    show_default=True,
    help="Keep solvent",
)
def show_pymol(filename, mask, mask2, sphere, stick, keep_solvent):
    """
    Use PYMOL to show H5/PDB/PARMED with mask

    \b
    Default: no solvent
    Optional two masks (green/purple). Atom mask selections, see:
      https://github.com/redesignScience/easytrajh5#atom-selection-language
    """
    easytrajh5.show.pymol(**locals())


@h5.command(no_args_is_help=True)
@click.argument("filename")
@click.argument("parmed")
@click.argument("mask", required=False)
@click.argument("mask2", required=False)
@click.option("-s", "--sphere", is_flag=True, help="Show all atoms as spheres")
@click.option("-t", "--stick", is_flag=True, help="Show all atoms as sticks")
@click.option(
    "-k",
    "--keep-solvent",
    is_flag=True,
    default=False,
    show_default=True,
    help="Keep solvent",
)
def show_chimera(filename, parmed, mask, mask2, sphere, stick, keep_solvent):
    """
    Use CHIMERA to show H5/PDB/PARMED with mask, needs PARMED

    \b
    Default: no solvent
    Optional two masks (green/purple). Atom mask selections, see:
      https://github.com/redesignScience/easytrajh5#atom-selection-language
    """
    easytrajh5.show.chimera(**locals())


@h5.command(no_args_is_help=True)
@click.argument("filename")
@click.argument("mask", required=False)
@click.argument("mask2", required=False)
@click.option("-s", "--sphere", is_flag=True, help="Show all atoms as spheres")
@click.option(
    "-k",
    "--keep-solvent",
    is_flag=True,
    default=False,
    show_default=True,
    help="Keep solvent",
)
def show_vmd(filename, mask, mask2, sphere, keep_solvent):
    """
    Use VMD to show H5/PDB/PARMED with mask

    \b
    Default: no solvent
    Optional two masks (green/purple). Atom mask selections, see:
      https://github.com/redesignScience/easytrajh5#atom-selection-language
    """
    easytrajh5.show.vmd(**locals())


if __name__ == "__main__":
    h5()
