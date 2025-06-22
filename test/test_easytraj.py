import filecmp
import importlib.util

from mdtraj import load_hdf5
from path import Path
import numpy

from easytrajh5.select import get_n_residue_of_mask
from easytrajh5.traj import EasyTrajH5File

this_dir = Path(__file__).absolute().parent
parmed = this_dir / "aaa-dry.parmed"
h5 = this_dir / "aaa_trajectory.h5"


def test_file_insertion():
    with EasyTrajH5File(h5) as f:
        f.insert_file_to_dataset("parmed", parmed)
    new_parmed = str(parmed).replace(".parmed", ".new.parmed")
    with EasyTrajH5File(h5) as f:
        f.extract_file_from_dataset("parmed", new_parmed)
    assert filecmp.cmp(parmed, new_parmed)


def test_get_parmed():
    h5 = this_dir / "aaa_trajectory.h5"
    pmd = EasyTrajH5File(h5).get_topology_parmed()
    n_solvent = get_n_residue_of_mask(pmd, "solvent")
    assert n_solvent > 0

    pmd = EasyTrajH5File(h5, atom_mask="not {solvent}").get_topology_parmed()
    n_solvent = get_n_residue_of_mask(pmd, "solvent")
    assert n_solvent == 0

    pmd = EasyTrajH5File(h5, atom_mask="mdtraj name CA").get_topology_parmed()
    assert len(pmd.atoms) == 3


def test_get_frames():
    h5 = this_dir / "aaa_trajectory.h5"
    f = EasyTrajH5File(h5)
    traj = f.read_as_traj()
    i_frames = [0, 5, 10]
    frames = f.read_frame_slice_as_traj(i_frames, None)
    for i_slice, i_frame in enumerate(i_frames):
        assert numpy.array_equal(traj.xyz[i_frame], frames.xyz[i_slice])


def test_blobs():
    with open(parmed, "rb") as f:
        blob = f.read()
    with EasyTrajH5File(h5) as f:
        f.set_bytes_dataset("parmed", blob)
    with EasyTrajH5File(h5) as f:
        reread_blob = f.get_bytes_dataset("parmed")
    assert blob == reread_blob


def test_json():
    stuff = {"haha": "ccc"}
    h5_copy = this_dir / "aaa_trajectory_j.h5"
    h5.copy(this_dir / "aaa_trajectory_j.h5")
    with EasyTrajH5File(h5_copy) as f:
        f.set_json_dataset("json_stuff", stuff)
    with EasyTrajH5File(h5_copy) as f:
        reread_stuff = f.get_json_dataset("json_stuff")
    assert stuff == reread_stuff


def test_easyh5_compatible_with_mdtraj():
    pytables_spec = importlib.util.find_spec("pytables")
    if pytables_spec is not None:
        mdtraj_traj = load_hdf5(h5)
        easy_traj = EasyTrajH5File(h5).read_as_traj()
        assert mdtraj_traj == easy_traj


def test_read_write_easyh5():
    easy_h5 = EasyTrajH5File(this_dir / "aaa_trajectory.h5")
    easy_traj = easy_h5.read_as_traj()

    keys = easy_h5.get_dataset_keys()
    traj_keys = [
        "coordinates",
        "time",
        "cell_lengths",
        "cell_angles",
        "velocities",
        "kineticEnergy",
        "potentialEnergy",
        "temperature",
        "alchemicalLambda",
    ]
    keys = [k for k in keys if k in traj_keys]

    n = easy_h5.get_n_frame()

    write_h5 = EasyTrajH5File(this_dir / "aaa_trajectory_w.h5", "w")
    write_h5.topology = easy_h5.topology

    for i in range(n):
        kwargs = {}
        for key in keys:
            kwargs[key] = easy_h5.get_dataset(key)[i]
        write_h5.write(**kwargs)
    write_h5.close()

    easy_h5.close()

    reread_h5 = EasyTrajH5File(this_dir / "aaa_trajectory_w.h5")
    test_mdtraj = reread_h5.read_as_traj()

    assert test_mdtraj == easy_traj


def test_h5_select_mask():
    h5 = this_dir / "aaa_trajectory.h5"
    with EasyTrajH5File(h5) as f:
        i_atoms = f.select_mask("amber @C")
        i_residues = f.select_mask_residues("amber @C")

        assert len(i_atoms) > 0
        assert numpy.array_equal(i_residues, [0, 1, 2])


