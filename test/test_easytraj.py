import filecmp
import shutil

import numpy
from mdtraj import Trajectory, load_hdf5
from path import Path
from pydash import py_

from easytrajh5.h5 import EasyH5
from easytrajh5.traj import EasyTrajH5File

THIS_DIR = Path(__file__).abspath().parent


def test_easyh5():
    # test data, 3 x frames of 4 x atoms of 3 coordinates
    conformations = numpy.array(
        [
            [
                [-0.12959832, -0.17306298, 1.05702889],
                [-0.16209413, -0.08169563, 1.16963661],
                [-0.08752456, -0.07581843, 1.28565717],
                [-0.13932393, -0.0190703, 1.39595604],
            ],
            [
                [-0.08647861, -0.21649601, 1.01632965],
                [-0.10281074, -0.11422855, 1.11324763],
                [-0.04765439, -0.12429251, 1.23844564],
                [-0.09414303, -0.0367968, 1.33750021],
            ],
            [
                [-0.06129837, -0.17449807, 1.04167926],
                [-0.06665925, -0.05670742, 1.12470198],
                [-0.02210106, -0.07005127, 1.2537868],
                [-0.0995131, -0.01841196, 1.35774732],
            ],
        ]
    )

    frame_shape = conformations.shape[1:]
    print(conformations.shape, frame_shape, conformations.dtype)
    stuff = {"haha": "ccc"}
    i_ligand_atoms = list(range(10))
    conf_h5 = THIS_DIR / "conformations.h5"

    with EasyH5(conf_h5, "w") as f:
        f.set_array_dataset("i_ligand_atoms", i_ligand_atoms)
        f.create_extendable_dataset("conformation", frame_shape, numpy.float32)
        f.create_extendable_dataset("another/conformation", frame_shape, numpy.float32)
        f.set_json_dataset("stuff", stuff)

    for c in conformations:
        with EasyH5(conf_h5, "a") as f:
            f.extend_dataset("conformation", [c])

    with EasyH5(conf_h5, mode="a") as f:
        assert f.get_dataset("conformation").shape == conformations.shape
        assert numpy.array_equal(f.get_dataset("i_ligand_atoms")[:], i_ligand_atoms)
        assert f.get_json_dataset("stuff") == stuff

    with EasyH5(THIS_DIR / "aaa_trajectory.h5", mode="r") as f:
        assert f.has_dataset("coordinates")
        assert f.get_dataset("coordinates").shape[2] == 3


def test_file_insertion():
    parmed = py_.head(THIS_DIR.glob("*.parmed"))
    traj = EasyTrajH5File(THIS_DIR / "aaa_trajectory.h5")
    traj.insert_file_to_dataset("parmed", parmed)
    traj.close()

    new_parmed = str(parmed).replace(".parmed", ".new.parmed")
    new_traj = EasyTrajH5File(THIS_DIR / "aaa_trajectory.h5")
    new_traj.extract_file_from_dataset("parmed", new_parmed)
    new_traj.close()

    assert filecmp.cmp(parmed, new_parmed)


def test_blobs():
    parmed = THIS_DIR.glob("*.parmed")[0]
    with open(parmed, "rb") as f:
        blob = f.read()

    traj = EasyTrajH5File(THIS_DIR / "aaa_trajectory.h5")
    traj.set_bytes_dataset("parmed", blob)
    traj.close()

    new_traj = EasyTrajH5File(THIS_DIR / "aaa_trajectory.h5")
    reread_blob = new_traj.get_bytes_dataset("parmed")
    new_traj.close()

    assert blob == reread_blob


def test_json():
    stuff = {"haha": "ccc"}
    shutil.copy(THIS_DIR / "aaa_trajectory.h5", THIS_DIR / "aaa_trajectory_j.h5")
    traj = EasyTrajH5File(THIS_DIR / "aaa_trajectory_j.h5")
    traj.set_json_dataset("json_stuff", stuff)
    traj.close()

    new_traj = EasyTrajH5File(THIS_DIR / "aaa_trajectory_j.h5")
    reread_stuff = new_traj.get_json_dataset("json_stuff")
    new_traj.close()

    assert stuff == reread_stuff


def test_easyh5_compatible_with_mdtraj():
    mdtraj_traj: Trajectory = load_hdf5(THIS_DIR / "aaa_trajectory.h5")

    easy_h5 = EasyTrajH5File(THIS_DIR / "aaa_trajectory.h5")
    easy_traj = easy_h5.read_as_traj()

    assert mdtraj_traj == easy_traj


def test_read_write_easyh5():
    easy_h5 = EasyTrajH5File(THIS_DIR / "aaa_trajectory.h5")
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

    write_h5 = EasyTrajH5File(THIS_DIR / "aaa_trajectory_w.h5", "w")
    write_h5.topology = easy_h5.topology

    for i in range(n):
        kwargs = {}
        for key in keys:
            kwargs[key] = easy_h5.get_dataset(key)[i]
        write_h5.write(**kwargs)
    write_h5.close()

    easy_h5.close()

    reread_h5 = EasyTrajH5File(THIS_DIR / "aaa_trajectory_w.h5")
    test_mdtraj = reread_h5.read_as_traj()

    assert test_mdtraj == easy_traj
