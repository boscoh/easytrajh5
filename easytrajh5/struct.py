import logging
import operator
import pickle
from pathlib import Path

import mdtraj
import parmed
from mdtraj.core import element as elem
from parmed import unit

from .fs import tic, toc
from .pdb import remove_model_lines
from .select import select_mask

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


def get_parmed(parmed_or_pdb: str) -> parmed.Structure:
    """
    :param parmed_or_pdb: str - either .parmed or .pdb
    """
    parmed_or_pdb = str(parmed_or_pdb)
    if Path(parmed_or_pdb).suffix == ".parmed":
        return load_parmed(parmed_or_pdb)
    elif Path(parmed_or_pdb).suffix == ".pdb":
        # Check for issue where mdtraj saves MODEL 0, which throws error in parmed
        remove_model_lines(parmed_or_pdb)
        return parmed.load_file(parmed_or_pdb)
    else:
        raise ValueError(
            f"Can't process {parmed_or_pdb} of type {Path(parmed_or_pdb.suffix)}, only .pdb and .parmed"
        )


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


def slice_parmed(pmd: parmed.Structure, atom_mask: str) -> (parmed.Structure, [int]):
    logger.info(tic("parsing atom mask"))
    i_atoms = select_mask(pmd, atom_mask, is_fail_on_empty=False)
    logger.info(toc())

    # handle parmed bug!
    if len(i_atoms) == len(pmd.atoms):
        return pmd, i_atoms

    logger.info(tic("slicing parmed"))
    sliced_pmd = pmd[i_atoms]
    logger.info(toc())

    return sliced_pmd, i_atoms


def slice_mdtraj_topology(mdtraj_top: mdtraj.Topology, atom_mask: str) -> (mdtraj.Topology, [int]):
    logger.info(tic("building parmed"))
    pmd = get_parmed_from_openmm(mdtraj_top.to_openmm())
    logger.info(toc())

    sliced_pmd, i_atoms = slice_parmed(pmd, atom_mask)

    logger.info(tic("converting to mdtraj topology"))
    sliced_top = mdtraj.Topology.from_openmm(sliced_pmd.topology)
    logger.info(toc())

    return sliced_top, i_atoms


def get_dict_from_mdtraj_topology(topology):
    try:
        topology_dict = {"chains": [], "bonds": []}

        for chain in topology.chains:
            chain_dict = {"residues": [], "index": int(chain.index)}
            for residue in chain.residues:
                residue_dict = {
                    "index": int(residue.index),
                    "name": str(residue.name),
                    "atoms": [],
                    "resSeq": int(residue.resSeq),
                    "segmentID": str(residue.segment_id),
                }

                for atom in residue.atoms:
                    try:
                        element_symbol_string = str(atom.element.symbol)
                    except AttributeError:
                        element_symbol_string = ""

                    residue_dict["atoms"].append(
                        {
                            "index": int(atom.index),
                            "name": str(atom.name),
                            "element": element_symbol_string,
                        }
                    )
                chain_dict["residues"].append(residue_dict)
            topology_dict["chains"].append(chain_dict)

        for atom1, atom2 in topology.bonds:
            topology_dict["bonds"].append([int(atom1.index), int(atom2.index)])

        return topology_dict

    except AttributeError as e:
        raise AttributeError(
            "topology_object fails to implement the"
            "chains() -> residue() -> atoms() and bond() protocol. "
            "Specifically, we encountered the following %s" % e
        )


def get_mdtraj_topology_from_dict(topology_dict) -> mdtraj.Topology:
    topology = mdtraj.Topology()

    for chain_dict in sorted(topology_dict["chains"], key=operator.itemgetter("index")):
        chain = topology.add_chain()
        for residue_dict in sorted(
            chain_dict["residues"], key=operator.itemgetter("index")
        ):
            try:
                resSeq = residue_dict["resSeq"]
            except KeyError:
                resSeq = None
                logger.warning(
                    "No resSeq information found in HDF file, defaulting to zero-based indices"
                )
            try:
                segment_id = residue_dict["segmentID"]
            except KeyError:
                segment_id = ""
            residue = topology.add_residue(
                residue_dict["name"], chain, resSeq=resSeq, segment_id=segment_id
            )
            for atom_dict in sorted(
                residue_dict["atoms"], key=operator.itemgetter("index")
            ):
                try:
                    element = elem.get_by_symbol(atom_dict["element"])
                except KeyError:
                    element = elem.virtual
                topology.add_atom(atom_dict["name"], element, residue)

    atoms = list(topology.atoms)
    for index1, index2 in topology_dict["bonds"]:
        topology.add_bond(atoms[index1], atoms[index2])

    return topology


