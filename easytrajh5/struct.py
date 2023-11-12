import itertools
import logging
import pickle

import mdtraj
import parmed
import pydash as py_
from parmed import unit
from path import Path

from .pdb import remove_model_lines

logger = logging.getLogger(__name__)

__doc__ = """
Useful transforms for parmed.Structure, mdtraj.Trajectory, and OpenMM
"""


def dump_parmed(pmd: parmed.Structure, fname: str):
    with open(fname, "wb") as handle:
        pickle.dump(file=handle, obj=pmd.__getstate__())


def load_parmed(fname: str) -> parmed.Structure:
    with open(fname, "rb") as handle:
        pmd = parmed.structure.Structure()
        pmd.__setstate__(pickle.load(file=handle))
    return pmd


def get_parmed_from_pdb(pdb: str) -> parmed.Structure:
    """
    Reads pdb with some sanity checks for model lines

    :param pdb: str - either .parmed or .pdb
    """
    suffix = Path(pdb).ext.lower()
    if not suffix == ".pdb":
        raise ValueError(f"Can't process {pdb} of type {suffix}, only .pdb")
    # Check for issue where mdtraj saves MODEL 0, which throws error in parmed
    remove_model_lines(pdb)
    return parmed.load_file(pdb)


def get_parmed_from_parmed_or_pdb(pdb_or_parmed: str) -> parmed.Structure:
    """
    :param pdb_or_parmed: str - either .parmed or .pdb
    """
    suffix = Path(pdb_or_parmed).ext
    if suffix == ".pdb":
        pmd = get_parmed_from_pdb(pdb_or_parmed)
    elif suffix == ".parmed":
        pmd = load_parmed(pdb_or_parmed)
    else:
        raise ValueError(f"Can't process {pdb_or_parmed} of type {suffix}")
    return pmd


def get_parmed_from_mdtraj(traj: mdtraj.Trajectory, i_frame=0) -> parmed.Structure:
    return parmed.openmm.load_topology(traj.top.to_openmm(), xyz=traj.xyz[i_frame])


def get_parmed_from_openmm(openmm_topology, openmm_positions=None) -> parmed.Structure:
    return parmed.openmm.load_topology(openmm_topology, xyz=openmm_positions)


def get_mdtraj_from_parmed(pmd: parmed.Structure) -> mdtraj.Trajectory:
    return mdtraj.Trajectory(
        xyz=pmd.coordinates / 10, topology=mdtraj.Topology.from_openmm(pmd.topology)
    )


def get_mdtraj_from_openmm(openmm_topology, openmm_positions) -> mdtraj.Trajectory:
    if unit.is_quantity(openmm_positions):
        openmm_positions = openmm_positions.value_in_unit(unit.nanometer)
    mdtraj_topology = mdtraj.Topology.from_openmm(openmm_topology)
    return mdtraj.Trajectory(topology=mdtraj_topology, xyz=openmm_positions)


def calc_residue_contacts_with_mdtraj(
    traj, i_residues1, i_residues2, cutoff_nm=None, max_n_residue=None
) -> [int]:
    """
    :return: [int] - indices of closest residues to ligand
    """
    # Generate pairs of residue indices [[i_lig1, i_res1], [i_lig1, i_res2]....]
    pairs = list(itertools.product(i_residues1, i_residues2))

    # Calculate distances as nx1 numpy.array and pairs is nx2 numpy.array
    # periodic=False turns off period cell correction
    distances, pairs = mdtraj.compute_contacts(
        traj, contacts=pairs, scheme="closest-heavy", periodic=False
    )

    # Get sorted top_entries list of contact residues
    top_entries = [(d, pair[1]) for d, pair in zip(distances[0], pairs)]
    top_entries = py_.sort_by(top_entries, lambda e: e[0])
    if max_n_residue:
        top_entries = top_entries[:max_n_residue]
    if cutoff_nm:
        top_entries = py_.filter_(top_entries, lambda e: e[0] <= cutoff_nm)

    return [e[1] for e in top_entries]


