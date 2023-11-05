import mdtraj
import mdtraj as md
import numpy
from mdtraj.utils import ilen, in_units_of


def filter_for_atom_lines(lines):
    return [l for l in lines if l.startswith("ATOM") or l.startswith("HETATM")]


def get_pdb_lines(
    mdtraj_topology,
    positions_angstroms,
    unitcell_lengths=None,
    unitcell_angles=None,
    bfactors=None,
    model_index=None,
    is_write_header=True,
):
    result = []

    if is_write_header and unitcell_lengths is not None and unitcell_angles is not None:
        box = list(unitcell_lengths) + list(unitcell_angles)
        assert len(box) == 6
        result.append(
            "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1 " % tuple(box)
        )

    if ilen(mdtraj_topology.atoms) != len(positions_angstroms):
        raise ValueError("The number of positions must match the number of atoms")
    if numpy.any(numpy.isnan(positions_angstroms)):
        raise ValueError("Particle position is NaN")
    if numpy.any(numpy.isinf(positions_angstroms)):
        raise ValueError("Particle position is infinite")

    if bfactors is None:
        bfactors = ["{0:5.2f}".format(0.0)] * len(positions_angstroms)
    else:
        if (numpy.max(bfactors) >= 100) or (numpy.min(bfactors) <= -10):
            raise ValueError("bfactors must be in (-10, 100)")
        bfactors = ["{0:5.2f}".format(b) for b in bfactors]

    if model_index is not None:
        result.append("MODEL     %4d" % model_index)

    i_atom = 1
    i_pos = 0
    _chain_names = [chr(ord("A") + i) for i in range(26)]
    for chainIndex, chain in enumerate(mdtraj_topology.chains):
        chainName = _chain_names[chainIndex % len(_chain_names)]
        residues = list(chain.residues)
        for resIndex, res in enumerate(residues):
            if len(res.name) > 3:
                resName = res.name[:3]
            else:
                resName = res.name
            for atom in res.atoms:
                if (
                    len(atom.name) < 4
                    and atom.name[:1].isalpha()
                    and (atom.element is None or len(atom.element.symbol) < 2)
                ):
                    atomName = " " + atom.name
                elif len(atom.name) > 4:
                    atomName = atom.name[:4]
                else:
                    atomName = atom.name
                coords = positions_angstroms[i_pos]
                if atom.element is not None:
                    symbol = atom.element.symbol
                else:
                    symbol = " "
                if atom.serial is not None and len(mdtraj_topology._chains) < 2:
                    # We can't do this for more than 1 chain
                    # to prevent issue 1611
                    atomSerial = atom.serial
                else:
                    atomSerial = i_atom
                line = (
                    "ATOM  %5d %-4s %3s %1s%4d    %s%s%s  1.00 %5s      %-4s%2s  "
                    % (  # Right-justify atom symbol
                        atomSerial % 100000,
                        atomName,
                        resName,
                        chainName,
                        (res.resSeq) % 10000,
                        _format_83(coords[0]),
                        _format_83(coords[1]),
                        _format_83(coords[2]),
                        bfactors[i_pos],
                        atom.segment_id[:4],
                        symbol[-2:],
                    )
                )
                assert len(line) == 80, "Fixed width overflow detected"
                result.append(line)
                i_pos += 1
                i_atom += 1
            if resIndex == len(residues) - 1:
                result.append(
                    "TER   %5d      %3s %s%4d"
                    % (atomSerial + 1, resName, chainName, res.resSeq)
                )
                i_atom += 1

    if model_index is not None:
        result.append("ENDMDL")

    return result


def get_pdb_lines_of_traj_frame(
    traj: mdtraj.Trajectory,
    i_frame=0,
    bfactors=None,
    model_index=None,
    is_write_header=True,
):
    """
    Adapted from mdtraj.formats.pdb.pdbfile.py
    """
    positions = in_units_of(traj._xyz[i_frame], traj._distance_unit, "angstroms")

    if traj.unitcell_lengths is None and traj.unitcell_angles is None:
        unitcell_lengths = None
        unitcell_angles = None
    elif traj.unitcell_lengths is None:
        raise ValueError("unitcell_lengths is missing")
    elif traj.unitcell_angles is None:
        raise ValueError("unitcell_angles is missing")
    else:
        unitcell_lengths = traj.unitcell_lengths[i_frame]
        unitcell_lengths = in_units_of(
            unitcell_lengths, traj._distance_unit, "angstroms"
        )
        if not len(unitcell_lengths) == 3:
            raise ValueError("unitcell_lengths must be length 3")
        unitcell_angles = traj.unitcell_angles[i_frame]
        if not len(unitcell_angles) == 3:
            raise ValueError("unitcell_angles must be length 3")

    return get_pdb_lines(
        traj.topology,
        positions,
        unitcell_lengths,
        unitcell_angles,
        bfactors,
        model_index,
        is_write_header,
    )


def save_traj_frame_pdb(traj, pdb, i_frame=0, model_index=None):
    with open(pdb, "w") as f:
        for line in get_pdb_lines_of_traj_frame(
            traj, i_frame=i_frame, model_index=model_index
        ):
            f.write(line + "\n")


def _format_83(f):
    """Format a single float into a string of width 8, with ideally 3 decimal
    places of precision. If the number is a little too large, we can
    gracefully degrade the precision by lopping off some of the decimal
    places. If it's much too large, we throw a ValueError"""
    if -999.999 < f < 9999.999:
        return "%8.3f" % f
    if -9999999 < f < 99999999:
        return ("%8.3f" % f)[:8]
    raise ValueError(
        'coordinate "%s" could not be represnted ' "in a width-8 field" % f
    )


def remove_model_lines(pdb):
    """
    Remove the "MODEL" and "ENDMDL" lines so as not confuse parmed
    """
    has_model_line = False
    lines = []
    with open(pdb) as f:
        for l in f.readlines():
            if "MODEL" in l:
                has_model_line = True
            elif "ENDMDL" in l:
                break
            else:
                lines.append(l)
    if has_model_line:
        with open(pdb, "w") as f:
            f.write("".join(lines))


def _split_pdb_file(pdb_file, ligand_residue_name) -> (str, str):
    pdb_file = str(pdb_file)
    traj = md.load_pdb(pdb_file)

    protein_ixs = traj.top.select(f"not resname {ligand_residue_name}")
    ligand_ixs = traj.top.select(f"resname {ligand_residue_name}")

    if len(ligand_ixs) == 0:
        raise ValueError(
            f"No ligand with residue name {ligand_residue_name} found in {pdb_file}"
        )

    protein_pdb_file = pdb_file.replace(".pdb", "_protein.pdb")
    protein = traj.atom_slice(protein_ixs)
    protein.save(protein_pdb_file)

    ligand_pdb_file = pdb_file.replace(".pdb", f"_{ligand_residue_name}.pdb")
    traj = traj.atom_slice(ligand_ixs)
    traj.save(ligand_pdb_file)

    return protein_pdb_file, ligand_pdb_file
