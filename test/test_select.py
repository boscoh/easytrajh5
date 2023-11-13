from copy import deepcopy

import numpy
from path import Path

from pydash import py_
from easytrajh5.select import parse_number_list, select_residue_contacts
from easytrajh5.select import select_mask, is_integer_list
from easytrajh5.struct import get_mdtraj_from_parmed, get_parmed_from_pdb

this_dir = Path(__file__).parent
pdb = this_dir / "pdb" / "hif2a__1__dock.pdb"


def is_equal_unsorted_array(a, b):
    a = deepcopy(a)
    b = deepcopy(b)
    a.sort()
    b.sort()
    return numpy.array_equal(a, b)


def test_number_list():
    for num_s, num_list in [
        ("1 2 3", [1, 2, 3]),
        ("50-100", list(range(50, 101))),
        ("2-3,5-7", [2, 3, 5, 6, 7]),
        ("7-9,11,13-15 19,21", [7, 8, 9, 11, 13, 14, 15, 19, 21]),
    ]:
        parse_list = parse_number_list(num_s)
        assert is_equal_unsorted_array(parse_list, num_list)

        parse_list = parse_number_list(num_list)
        assert is_equal_unsorted_array(parse_list, num_list)


def test_protein_ligand_selections():
    pmd = get_parmed_from_pdb(pdb)

    l_resname = "UNL"
    ligand_mask = f"amber :{l_resname}" if l_resname else "ligand"
    complex_sel = select_mask(pmd, "merge {protein} {" + ligand_mask + "}")
    protein_sel = select_mask(pmd, "protein")
    ligand_sel = select_mask(pmd, ligand_mask)

    assert len(ligand_sel) != len(complex_sel)
    assert len(ligand_sel) + len(protein_sel) == len(complex_sel)


def test_mdtraj_and_parmed_indexing():
    pmd = get_parmed_from_pdb(pdb)
    traj = get_mdtraj_from_parmed(pmd)
    for pmd_res, traj_res in zip(pmd.residues, traj.top.residues):
        assert pmd_res.number == traj_res.resSeq
        assert pmd_res.idx == traj_res.index
        assert pmd_res.name == traj_res.name
    for pmd_atom, traj_atom in zip(pmd.atoms, traj.top.atoms):
        assert pmd_atom.idx == traj_atom.index


def test_amber_indexing():
    resi_list = [13, 15, 17]
    mask = "amber :" + ",".join(map(str, resi_list))

    pmd = get_parmed_from_pdb(pdb)
    i_atoms = select_mask(pmd, mask)
    found_resi_list = py_.uniq([pmd.atoms[i].residue.number for i in i_atoms])
    assert is_equal_unsorted_array(resi_list, found_resi_list)

    traj = get_mdtraj_from_parmed(pmd)
    atoms = list(traj.top.atoms)
    found_resi_list = py_.uniq([atoms[i].residue.resSeq for i in i_atoms])
    assert is_equal_unsorted_array(resi_list, found_resi_list)


def test_resi_selections():
    resi_list = [13, 15, 17]
    num_list = " ".join(map(str, resi_list))

    pmd = get_parmed_from_pdb(pdb)
    i_atoms = select_mask(pmd, f"resi {num_list}")

    found_resi_list = []
    for i in i_atoms:
        found_resi_list.append(pmd.atoms[i].residue.idx)
    found_resi_list = py_.uniq(found_resi_list)

    assert is_equal_unsorted_array(resi_list, found_resi_list)


def test_atom_selections():
    pmd = get_parmed_from_pdb(pdb)
    atom_list = [13, 15, 17]
    num_list = " ".join(map(str, atom_list))

    i_atoms = select_mask(pmd, f"atom {num_list}")
    assert is_equal_unsorted_array(i_atoms, atom_list)

    for i in i_atoms:
        assert pmd.atoms[i].idx == i

    i_atoms = select_mask(pmd, atom_list)
    assert is_equal_unsorted_array(i_atoms, atom_list)


def test_integer_list_checker():
    assert is_integer_list([numpy.int32(1), 2, 3])
    assert is_integer_list(numpy.array([1, 2, 3]))
    assert is_integer_list([1, 2, 3])
    assert not is_integer_list(1)
    assert not is_integer_list([1, 2, "3"])


def test_contact():
    pmd = get_parmed_from_pdb(pdb)

    i_contact_atoms = select_residue_contacts(pmd, "UNL")

    i_contact_by_mask_atoms = select_mask(pmd, "near UNL")
    assert is_equal_unsorted_array(i_contact_atoms, i_contact_by_mask_atoms)

    i_contact_by_mask_atoms = select_mask(pmd, "pocket")
    assert is_equal_unsorted_array(i_contact_atoms, i_contact_by_mask_atoms)


if __name__ == "__main__":
    test_amber_indexing()
