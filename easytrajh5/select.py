import itertools
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

import mdtraj
import parmed
import pydash as py_
from parmed.amber import AmberMask

from .fs import load_yaml

data_dir = Path(__file__).parent
resnames_by_keyword = load_yaml(data_dir / "data" / "select.yaml")


def get_resnames(keyword):
    return resnames_by_keyword[keyword]


def select_mask(pmd, mask, temp_pdb="temp.pdb", is_fail_on_empty=True):
    """
    Selects atom based on a selection string which is based on a combination
    atom selection langauge.

    There are several different types of selection modes strings, where
    the first word is oftend used a mode selector.

    - ommtk keywords
        "protein pocket ligand noh" - selects protein, pocket, ligand atoms and skips waters
    - ommtk resid selections
        "resid A 0 10" - selects atoms in the first and 11th residue in chain A
    - atom 0-indexed indices
        "atom 0 55 43" - selects the first, 56th, and 44th atom
    - AMBER-style atom selection https://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax
        "amber :ALA,LYS" - selects all alanine and lysine residues
    - MDTraj-style atom selection
        "mdtraj protein and water" - selects protein and water
    - furthermore, selections can be combined with set operators ("not", "intersect", "merge", "diff"),
        "intersect {not {amber :ALA}} {protein}"
        "diff {protein} {not {amber :ALA}}"
        "not {resname LEU}"
        "merge {near BSM 8} {amber :ALA}"

    :parmam parmed_structure: parmed.Structure
    :param mask: Union[str, List[int]] - described above
    :return: [int]
    """

    def get_i_atoms_of_ast(ast):
        has_list = py_.some(ast, lambda node: isinstance(node, list))
        if not has_list:
            expr = " ".join(ast)
            value = process_expr(pmd, expr, temp_pdb)
            return value

        operator = ast[0]

        if operator == "not":
            if len(ast) != 2:
                raise ValueError("not: must have only 1 following expr in {}")
            i_atoms1 = get_i_atoms_of_ast(ast[1])
            return select_not_atoms(pmd, i_atoms1)

        elif operator == "merge":
            result = []
            for node in ast[1:]:
                result.extend(get_i_atoms_of_ast(node))
            return list(set(result))

        elif operator == "intersect":
            result = set(get_i_atoms_of_ast(ast[1]))
            for node in ast[2:]:
                i_atom_set = set(get_i_atoms_of_ast(node))
                result = result.intersection(i_atom_set)
            return list(result)

        elif operator == "diff":
            if len(ast) != 3:
                raise ValueError("diff: must have 2 expr in {}")
            i_atoms1 = get_i_atoms_of_ast(ast[1])
            i_atoms2 = get_i_atoms_of_ast(ast[2])
            return diff_list(i_atoms1, i_atoms2)

        else:
            raise ValueError(f"Operator {operator} not in [not diff merge intersect]")

    i_atoms = get_i_atoms_of_ast(parse_ast(mask))

    i_atoms.sort()

    res_indices = []
    for i in i_atoms:
        res_indices.append(pmd.atoms[i].residue.idx)
    res_indices = py_.uniq(res_indices)
    logger.info(
        f'select_mask "{mask}" -> {len(i_atoms)} atoms, {len(res_indices)} residues'
    )

    if is_fail_on_empty and not len(i_atoms):
        raise ValueError("Selection produced no atoms")

    return i_atoms


def process_expr(pmd, expr, temp_pdb):
    if expr.startswith("amber "):
        amber_mask = AmberMask(pmd, expr[6:])
        i_atoms = [i_atom for i_atom, mask in enumerate(amber_mask.Selection()) if mask]
    elif expr.startswith("mdtraj "):
        mdtraj_top = mdtraj.Topology.from_openmm(pmd.topology)
        return mdtraj_top.select(expr[6:]).tolist()
    elif expr.startswith("atom"):
        i_atoms = [int(x) for x in expr[5:].split()]
    elif expr.startswith("resid "):
        i_atoms = original_select_atoms(pmd, resid_selection=expr[6:])
    else:
        i_atoms = original_select_atoms(pmd, keyword_selection=expr, temp_pdb=temp_pdb)
    return i_atoms


def parse_ast(mask):
    """
    Parses an expression bounded by curly brackets into an abstract syntax tree
    e.g. 'a {{b c} {d}} e}' -> ['a', [[['b', 'c'], ['d'], 'e']]
    """

    def push(obj, a_list, depth):
        while depth:
            a_list = a_list[-1]
            depth -= 1
        a_list.append(obj)

    def parse_parentheses(tokens):
        groups = []
        depth = 0

        try:
            for token in tokens:
                if token == "{":
                    push([], groups, depth)
                    depth += 1
                elif token == "}":
                    depth -= 1
                else:
                    push(token, groups, depth)
        except IndexError:
            raise ValueError("Parentheses mismatch")

        if depth > 0:
            raise ValueError("Parentheses mismatch")
        else:
            return groups

    mask = mask.replace("{", " { ").replace("}", " } ")
    tokens = py_.filter_(mask.split(" "))
    return parse_parentheses(tokens)


def select_resid(pmd, resid_selection):
    """
    Selects atoms belonging to residues specified by the following residue string
    eg "A 11 12 13 14" where the first is the chain identifier and the following
    residues are also assumed to be joined by 'or' operator such that the above selection
    will return all atoms in residues 11, 12, 13 and 14 on chain A.

    :param pmd: parmed.Structure
    :param reside_selection: str
    :return: [int]
    """
    final_list = []

    # split up the selection syntax
    tokens = resid_selection.split()
    chain_selection = tokens[0]
    resid_list = [int(i) for i in tokens[1:]]

    # TODO check the order of chains in multichain proteins
    # TODO Currently assuming alphabet designation matches parmed order
    protein_chains = [
        i
        for i in pmd.topology.chains()
        if next(i.residues()).name in get_resnames("protein")
    ]
    for k, chain in enumerate(protein_chains):
        # convert letter rep of chain to ordinal number rep
        if ord(chain_selection.lower()) - 96 == k + 1:
            residues = [i for i in chain.residues()]
            res0 = residues[0]
            start_index = res0.index
            selected_residues = [
                i for i in residues if i.index - start_index in resid_list
            ]

            if len(selected_residues) == 0:
                raise ValueError(
                    "Could not find one of the residues {} on protein chain.".format(
                        resid_list
                    )
                )

            for i in selected_residues:
                for k in i.atoms():
                    final_list.append(k.index)
    return final_list


def calc_contacts(pdb, lig_resname, cutoff_nm=None, max_n_residue=6):
    """
    Finds the residues closests to ligand

    :param pdb: str
    :param lig_resname: str - name of ligand residue, assumes only one
    :param max_n_residue: int - maximum number of residues to use
    :return: [int] - indices of closest residues to ligand
    """
    traj: mdtraj.Trajectory = mdtraj.load_pdb(pdb)
    residues = list(traj.topology.residues)

    logger.info(f"calc_contacts {pdb} {lig_resname}")

    # Generate pairs of residue indices of ligand and protein residues
    i_ligands = [r.index for r in residues if r.name == lig_resname]
    if len(i_ligands) == 0:
        raise ValueError(f"no residue with resname={lig_resname}")
    i_residues = [r.index for r in residues if r.name in get_resnames("protein")]

    result = traj_calc_residue_contacts(
        traj, i_ligands, i_residues, cutoff_nm=cutoff_nm, max_n_residue=max_n_residue
    )
    print(f"closest {len(result)} residues")
    return result


def select_contact_residues(
        pmd, lig_resname="LIG", max_n_residue=6, temp_pdb="temp.pdb"
):
    """
    Finds the atoms of the residues closest to a ligand residue

    :return: [int]
    """
    if Path(temp_pdb).exists():
        Path(temp_pdb).unlink()
    pmd.save(temp_pdb)
    i_contact_residues = calc_contacts(temp_pdb, lig_resname, max_n_residue)
    resid_selection = "A " + " ".join([str(i) for i in i_contact_residues])
    return select_resid(pmd, resid_selection)


def select_resnames(pmd, resnames):
    return [a.idx for a in pmd.atoms if a.residue.name in resnames]


def select_element(pmd, element):
    return [a.idx for a in pmd.atoms if element in a.element_name]


def diff_list(list1, list2):
    return sorted(set(list1) - set(list2))


def select_not_atoms(parmed_structure, i_atoms):
    i_all_atoms = [a.idx for a in parmed_structure.atoms]
    return diff_list(i_all_atoms, i_atoms)


def original_select_atoms(
        pmd: parmed.Structure,
        keyword_selection=None,
        ligand_resname=None,
        resid_selection=None,
        temp_pdb=None,
) -> [int]:
    """
    Select a list of atom  that can be used to slice parmed objects, set restraints, or specify CVs.
    If both keyword and resid selection are applied both are returned (i.e. 'or' combining is assumed).

    :param keyword_selection:
        Accepts (in any order): 'ligand', 'protein', 'water', 'lipid'
        If more than one keyword is specified, it is assumed they are joined with "or"
        operation (i.e. 'ligand protein' will return both ligand and protein atom indices).
        Modifiers 'noh' and 'not' where 'noh' will exclude hydrogens and 'not' will invert the selection.
        The keyword 'ligand' will find the residue 'LIG', 'UNL', or whatever is in 'ligand_resname'
        The keyword 'pocket' will find the closest 6 residues to the ligand.
        The keyword 'near' will require a following resname, with an optional integer, e.g.:
            'near ATP'
            'near ATP 5'
        The keyword 'resname' will require a resname:
            'resname LEU'
        The keyword 'element' will require a resname:
            'element H'
    :param resid_selection: eg "A 11 12 13 14" where the first is the chain identifier and the following
        residues are also assumed to be joined by 'or' operator such that the above selection will return
        all atoms in residues 11, 12, 13 and 14 on chain A.

        If 'not' or 'noh' are selected the modifier will be applied to all selections (keyword or resid). As such
        >>> sel = original_select_atoms(pmd, keyword_selection='noh', resid_selection='A 11')
        will return non-hydrogen atoms in residue 11 on chain A.
    :param ligand_resname: by default 'LIG' and 'UNL' are assumed to be ligands resnames. If your ligand
        has a different residue name you can specify it with ligand_resname. Otherwise Amber style
        residue naming is assumed.
    """
    if keyword_selection == None and resid_selection == None:
        raise ValueError("Must specify either keyword selection or resid_selection")
    if ligand_resname:
        resnames_by_keyword["ligand"].append(
            ligand_resname
        )  # JM - added since the dict was not updating new ligand resnames
    result = []
    if resid_selection:
        result.extend(select_resid(pmd, resid_selection))
    if keyword_selection:
        tokens = [t for t in keyword_selection.split(" ") if t]
        is_noh = "noh" in tokens
        is_nosolvent = "nosolvent" in tokens
        allowed_keywords = list(resnames_by_keyword.keys()) + [
            "pocket",
            "near",
            "resname",
            "noh",
            "nosolvent",
        ]
        while len(tokens) > 0:
            keyword = tokens.pop(0)
            if keyword not in allowed_keywords:
                raise ValueError(f"Keyword {keyword} not in {allowed_keywords}")
            i_atoms = []
            if keyword in resnames_by_keyword:
                i_atoms = select_resnames(pmd, resnames_by_keyword[keyword])
            elif keyword == "pocket":
                i_atoms = select_contact_residues(pmd, "LIG", temp_pdb=temp_pdb)
            elif keyword == "near":
                if len(tokens) < 1:
                    raise ValueError("keyword near requires a resname argument")
                lig_resname = tokens.pop(0)
                if len(tokens) and tokens[0].isdigit():
                    n = int(tokens.pop(0))
                    i_atoms = select_contact_residues(
                        pmd,
                        lig_resname=lig_resname,
                        max_n_residue=n,
                        temp_pdb=temp_pdb,
                    )
                else:
                    i_atoms = select_contact_residues(
                        pmd, lig_resname=lig_resname, temp_pdb=temp_pdb
                    )
            elif keyword == "resname":
                if len(tokens) < 1:
                    raise ValueError("keyword res requires an argument")
                resname = tokens.pop(0)
                i_atoms = select_resnames(pmd, [resname])
                if not len(i_atoms):
                    logger.warning(
                        f"Warning: no atoms were found for resname={resname}"
                    )
            elif keyword not in ["noh", "nosolvent"]:
                raise ValueError(f"Can't parse '{keyword}'")
            result.extend(i_atoms)
        if is_noh:
            result = diff_list(result, select_element(pmd, "H"))
        if is_nosolvent:
            solvent = select_resnames(pmd, resnames_by_keyword["solvent"])
            result = diff_list(result, solvent)
    return sorted(set(result))


def traj_calc_residue_contacts(
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
