import importlib.util
import itertools
import logging
import pickle
import tempfile

import mdtraj
import parmed
import pydash as py_
from packaging import version
from parmed import unit
from parmed.topologyobjects import TrackedList
from path import Path

from .pdb import remove_model_lines

logger = logging.getLogger(__name__)

__doc__ = """
Useful transforms for parmed.Structure, mdtraj.Trajectory, and OpenMM.
Loads structures using convenient transforms and unit conversions.

We have found that saving parmed.Structure.__getstate__() as a pickle,
is a very efficient and fast way of saving cross-platform MD prep files.

For units:
  - mdtraj and openmm use nanometers
  - pdb & parmed use angstroms

"""


def patch_pmd(pmd):
    # Compatibility for parmed 3.4.3 => 4.*
    for a in pmd.atoms:
        for key in ["formal_charge", "hybridization", "aromatic"]:
            if not hasattr(a, key):
                setattr(a, key, None)
    return pmd


def dump_parmed(pmd: parmed.Structure, fname: str):
    pmd = patch_pmd(pmd)
    state = pmd.__getstate__()
    state["parmed_version"] = parmed.__version__
    if hasattr(pmd, "extra"):
        state["extra"] = pmd.extra
    with open(fname, "wb") as handle:
        pickle.dump(file=handle, obj=state)


def patch_parmed_state(state):
    # Compatibility for parmed 3.4.3 => 4.*
    if version.parse(parmed.__version__) > version.parse("3.4.3"):
        if "links" not in state:
            state["links"] = TrackedList()
    return state

def load_parmed(fname: str) -> parmed.Structure:
    with open(fname, "rb") as handle:
        state = pickle.load(file=handle)
    pmd = parmed.structure.Structure()
    pmd.__setstate__(patch_parmed_state(state))
    pmd = patch_pmd(pmd)
    if "extra" in state:
        pmd.extra = state["extra"]
    return pmd


def get_parmed_from_pdb(pdb: str) -> parmed.Structure:
    """
    Reads pdb with a sanity check for model lines that confuses parmed
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


def get_parmed_from_mdtraj_topology_and_positions(topology, positions_angstroms):
    openmm_spec = importlib.util.find_spec("openmm")
    if openmm_spec is not None:
        return parmed.openmm.load_topology(
            topology.to_openmm(), xyz=positions_angstroms
        )
    else:
        logger.info("no openmm: save temp.pdb for mdtraj->parmed")
        positions_nm = positions_angstroms / 10
        traj = mdtraj.Trajectory(positions_nm, topology)
        s = next(tempfile._get_candidate_names())
        temp_pdb = f"temp.{s}.pdb"
        traj.save_pdb(temp_pdb, force_overwrite=True)
        return get_parmed_from_pdb(temp_pdb)


def get_mdtraj_topology_from_pmd(pmd):
    openmm_spec = importlib.util.find_spec("openmm")
    if openmm_spec is not None:
        return mdtraj.Topology.from_openmm(pmd.topology)
    else:
        s = next(tempfile._get_candidate_names())
        temp_pdb = f"temp.{s}.pdb"
        pmd.save(temp_pdb, overwrite=True)
        return mdtraj.load_topology(temp_pdb)


def get_parmed_from_mdtraj(traj: mdtraj.Trajectory, i_frame=0) -> parmed.Structure:
    positions_angs = traj.xyz[i_frame] * 10
    return get_parmed_from_mdtraj_topology_and_positions(traj.top, positions_angs)


def get_parmed_from_openmm(
    openmm_topology, positions_angstroms=None
) -> parmed.Structure:
    """
    :param positions_angstroms: unit.Quantity(dist) | [float] in angstroms
    """
    return parmed.openmm.load_topology(openmm_topology, xyz=positions_angstroms)


def get_mdtraj_from_parmed(pmd: parmed.Structure) -> mdtraj.Trajectory:
    return mdtraj.Trajectory(
        xyz=pmd.coordinates / 10, topology=get_mdtraj_topology_from_pmd(pmd)
    )


def get_mdtraj_from_openmm(openmm_topology, positions_nm) -> mdtraj.Trajectory:
    if unit.is_quantity(positions_nm):
        positions_nm = positions_nm.value_in_unit(unit.nanometer)
    mdtraj_topology = mdtraj.Topology.from_openmm(openmm_topology)
    return mdtraj.Trajectory(topology=mdtraj_topology, xyz=positions_nm)


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


def slice_parmed(pmd: parmed.Structure, i_atoms: [int]) -> parmed.Structure:
    # This function avoids the issue where parmed expects a bit
    # mask for selections for full selections but atom indices otherwise
    if len(i_atoms) == len(pmd.atoms):
        return pmd
    result = pmd[i_atoms]
    if hasattr(pmd, "extra"):
        result.extra = pmd.extra
    return result
