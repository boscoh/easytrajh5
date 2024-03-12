import logging
import operator
import pickle
from collections import namedtuple
from typing import TypeVar, Sequence

import numpy
import parmed
from addict import Dict
from mdtraj import Topology, Trajectory
from mdtraj.core import element as elem
from mdtraj.reporters.basereporter import _BaseReporter
from mdtraj.utils import ensure_type, in_units_of
from mdtraj.utils.unitcell import lengths_and_angles_to_box_vectors
from path import Path
from pydash import py_

from easytrajh5.fs import load_yaml, dump_yaml
from easytrajh5.pdb import remove_model_lines
from easytrajh5.struct import patch_parmed_state, patch_pmd
from .fs import tic, toc
from .h5 import EasyH5File
from .select import select_mask
from .struct import (
    slice_parmed,
    get_parmed_from_mdtraj_topology_and_positions,
    get_mdtraj_topology_from_pmd,
)

logger = logging.getLogger(__name__)

_Slice = TypeVar("_Slice", bound=Sequence)


def check_is_int_list(a):
    if not isinstance(a, list):
        return False
    return all(isinstance(i, int) for i in a)


class EasyTrajH5File(EasyH5File):
    """
    EasyTrajH5File is an object used to conviently process mdtraj-type
    H5 files efficiently, using the internal indexing features
    of H5 via the h5py library.

    In particular, we provide a convenient atom selection language
    to pre-select for atoms before extracting the trajectory. This
    is used to initialize the object via the topology, and subsequent
    read methods will read only the selected atoms.

    As well, a special dry_cache option has been implemented that
    will store the dry topology ("not solvent") in a special cache
    inside the H5 file.  Accessing the reduced dry topology will be
    fast, as the H5 can skip the loading of the full topology.
    Subsequent frame reads will only reference the dry atoms.

    The trajectory can be returned as:

        - read_as_traj() - the whole trajectory
        - read_frame_as_traj() - a traj with one frame
        - read_frame_slice_as_traj() - a selection of frames

    As well, the API of mdtraj.formats.H5TrajectoryFile has been replicated,
    and EasyTrajH5File can serve as a backend of openmm._BaseReporter:

        - __init__(h5, mode)
        - @property.setter topology
        - distance_unit
        - flush (optional)
        - close
        - read
        - write
    """

    # Compatibility with mdtraj.formats.hdf5.H5TrajectoryFile.read
    # nanometers
    distance_unit = Trajectory._distance_unit

    fields = [
        Dict(
            key="coordinates",
            shape=(None, None, 3),
            units=distance_unit,
            traj_key="xyz",
        ),
        Dict(key="time", shape=(None,), units="picoseconds", traj_key="time"),
        Dict(
            key="cell_lengths",
            shape=(None, 3),
            units=distance_unit,
            traj_key="unitcell_lengths",
        ),
        Dict(
            key="cell_angles",
            shape=(None, 3),
            units="degrees",
            traj_key="unitcell_angles",
        ),
        Dict(key="velocities", shape=(None, None, 3), units="nanometers/picosecond"),
        Dict(key="kineticEnergy", shape=(None,), units="kilojoules_per_mole"),
        Dict(key="potentialEnergy", shape=(None,), units="kilojoules_per_mole"),
        Dict(key="temperature", shape=(None,), units="kelvin"),
        # reversed the mdtraj choice to rename this `lambda` in the h5
        Dict(key="alchemicalLambda", shape=(None,), units="dimensionless"),
    ]

    default_attrs = {
        "conventions": "Pande",
        "conventionVersion": "1.1",
        "program": "rseed",
        "programVersion": "0.1",
        "application": "rseed",
    }

    def __init__(
        self,
        fname: str,
        mode: str = "a",
        atom_mask: str = "",
        is_dry_cache: bool = False,
    ):
        logger.info(
            f"{self.__class__.__name__}: {fname=} {mode=} {atom_mask=} {is_dry_cache=}"
        )
        logger.info(tic("open connection"))
        super().__init__(fname, mode)
        logger.info(toc())

        # to reach here, we've successfully loaded an .h5 file
        if mode == "w":
            self.has_header = False
        elif mode in ["a", "r"]:
            self.has_header = self.has_dataset("coordinates")
        else:
            raise ValueError("mode must be one of ['r', 'w', 'a']")

        # stores the full topology
        self._topology = None
        self.topology_by_atom_hash = {}

        self.atom_mask = atom_mask
        self.atom_indices = None
        self.last_i_frame = None
        self.last_frame = None
        self.last_atom_indices = None

        if is_dry_cache:
            if atom_mask and atom_mask != "not {solvent}":
                raise ValueError(
                    "Can't set is_dry_cache=True and atom_mask at the same time"
                )
            atom_mask = "not {solvent}"
            logger.info(f"{mode=} {atom_mask=}")
            if mode == "r":
                logger.info("Warning: can't use dry_atom_cache in read-only mode")
            else:
                dry_atom_indices = self.load_dry_topology_from_cache()
                self.atom_indices = dry_atom_indices
                return

        if atom_mask:
            self.set_atom_mask(atom_mask)

    def set_atom_mask(self, atom_mask):
        logger.info("setting atom_indices from atom_mask")
        sliced_top, atom_indices = self.slice_topology(atom_mask)
        self.atom_indices = atom_indices
        atom_hash = tuple(atom_indices)
        self.topology_by_atom_hash[atom_hash] = sliced_top

    def fetch_topology(self, atom_indices=None) -> Topology:
        if atom_indices is not None:
            atom_hash = tuple(atom_indices)
            if atom_hash in self.topology_by_atom_hash:
                return self.topology_by_atom_hash[atom_hash]

        # Need to check if self._topology exists
        if self._topology is None:
            if not self.has_dataset("topology"):
                raise ValueError(f"No topology saved in {self.fname}")
            logger.info(tic("loading topology"))
            topology_dict = self.get_json_dataset("topology")
            self._topology = get_mdtraj_topology_from_dict(topology_dict)
            logger.info(toc())

        if atom_indices is not None:
            atom_hash = tuple(atom_indices)
            logger.info(tic(f"slicing topology {len(atom_indices)}"))
            topology = self._topology.subset(atom_indices)
            logger.info(toc())
            self.topology_by_atom_hash[atom_hash] = topology
            return topology
        else:
            return self._topology

    def load_dry_topology_from_cache(self):
        if self.has_dataset("dry_atoms"):
            logger.info(tic("reading dry_atoms"))
            dry_atom_indices = self.get_dataset("dry_atoms")[:]
            logger.info(toc())

            if len(dry_atom_indices) == 0:
                logger.info("dry_atoms=[] -> no dry_topology")
                dry_atom_indices = None
            elif self.has_dataset("dry_topology"):
                logger.info(tic("reading dry_topology"))
                topology_dict = self.get_json_dataset("dry_topology")
                dry_top = get_mdtraj_topology_from_dict(topology_dict)
                logger.info(toc())
                atom_hash = tuple(dry_atom_indices)
                self.topology_by_atom_hash[atom_hash] = dry_top
        else:
            # haven't tried calculating dry_topology yet
            dry_top, dry_atom_indices = self.slice_topology("not {solvent}")
            n_atom = sum(1 for _ in self.topology.atoms)
            if len(dry_atom_indices) == n_atom:
                logger.info(tic("no solvent => save dry_atoms=[]"))
                self.set_array_dataset("dry_atoms", numpy.array([]))
                logger.info(toc())
                dry_atom_indices = None
            else:
                logger.info(tic("saving dry_topology and dry_atoms"))
                self.set_json_dataset(
                    "dry_topology", get_dict_from_mdtraj_topology(dry_top)
                )
                self.set_array_dataset("dry_atoms", numpy.array(dry_atom_indices))
                self.flush()
                logger.info(toc())
                atom_hash = tuple(dry_atom_indices)
                self.topology_by_atom_hash[atom_hash] = dry_top

        return dry_atom_indices

    def get_topology_parmed(self, i_frame=0):
        positions_angstroms = self.get_dataset("coordinates")[i_frame] * 10
        if self.atom_indices is None:
            mdtraj_topology = self.topology
        else:
            logger.info("Select atom_mask")
            mdtraj_topology = self.fetch_topology(self.atom_indices)
            positions_angstroms = positions_angstroms[self.atom_indices]
        return get_parmed_from_mdtraj_topology_and_positions(
            mdtraj_topology, positions_angstroms
        )

    def get_full_topology_parmed(self, i_frame=0):
        positions_angstroms = self.get_dataset("coordinates")[i_frame] * 10
        mdtraj_topology = self.topology
        return get_parmed_from_mdtraj_topology_and_positions(
            mdtraj_topology, positions_angstroms
        )

    def get_parmed_from_dataset(self, i_frame=None):
        if not self.has_dataset("parmed"):
            return None

        pmd = parmed.Structure()
        state = pickle.loads(self.get_bytes_dataset("parmed"))
        pmd.__setstate__(patch_parmed_state(state))
        pmd = patch_pmd(pmd)

        if i_frame is not None:
            coordinates = self.get_dataset("coordinates")
            if i_frame < 0:
                i_frame = coordinates.shape[0] + i_frame
            pmd.positions = coordinates[i_frame] * 10.0
            box_vectors = self.read_parmed_box_vectors(i_frame)
            if box_vectors is not None:
                pmd.box_vectors = box_vectors

        return pmd

    def slice_topology(self, atom_mask: str) -> (Topology, [int]):
        logger.info(tic("building parmed"))
        positions_angstroms = self.get_dataset("coordinates")[0] * 10
        pmd = get_parmed_from_mdtraj_topology_and_positions(
            self.topology, positions_angstroms
        )
        logger.info(toc())

        logger.info(tic("select mask"))
        i_atoms = select_mask(pmd, atom_mask, is_fail_on_empty=False)
        logger.info(toc())

        logger.info(tic("slicing parmed"))
        sliced_pmd = slice_parmed(pmd, i_atoms)
        logger.info(toc())

        logger.info(tic("converting to mdtraj topology"))
        sliced_top = get_mdtraj_topology_from_pmd(sliced_pmd)
        logger.info(toc())

        return sliced_top, i_atoms

    def select_mask(self, mask):
        pmd = self.get_full_topology_parmed()
        return select_mask(pmd, mask, is_fail_on_empty=False)

    def select_mask_residues(self, mask):
        pmd = self.get_full_topology_parmed()
        i_atoms = select_mask(pmd, mask, is_fail_on_empty=False)
        i_residues = [pmd[i_atom].residue.idx for i_atom in i_atoms]
        return py_.sort(py_.uniq(i_residues))

    @property
    def topology(self) -> Topology:
        if self._topology is None:
            self.fetch_topology()
        return self._topology

    @topology.setter
    def topology(self, topology):
        """Compatibility with mdtraj.formats.hdf5.H5TrajectoryFile.topology="""
        self._topology = topology
        self.set_json_dataset("topology", get_dict_from_mdtraj_topology(topology))
        self.handle.flush()

    def get_n_frame(self):
        return self.get_dataset("coordinates").shape[0]

    def get_n_atom(self):
        return self.get_dataset("coordinates").shape[1]

    def __len__(self):
        return self.get_n_frame()

    def read_atom_dataset_progressively(
        self, key, frame_slice, atom_indices=slice(None)
    ) -> numpy.ndarray:
        """Returns numpy.ndarrary[float] of [n_frame x <shape of frame of dataset>]"""

        stride = frame_slice.step or 1
        start = frame_slice.start or 0

        dataset = self.get_dataset(key)
        n_frame_total = frame_slice.stop - start
        n_atom = (
            len(atom_indices) if isinstance(atom_indices, list) else dataset.shape[1]
        )
        n_frame_in_chunk = dataset.chunks[0]
        n_frame = n_frame_total // stride
        offset = 1 if n_frame % stride != 0 else 0
        result = numpy.empty((n_frame + offset, n_atom, 3), dataset.dtype)
        # Try this n_frame_in_page first, and let errors force us to reduce it
        n_frame_of_page = n_frame_total

        i_frame_of_page = start
        i_frame_in_result = 0
        while True:
            try:
                frame_slice = slice(
                    i_frame_of_page, i_frame_of_page + n_frame_of_page, stride
                )
                data = dataset[frame_slice, atom_indices, :]
            except IOError as ioe:
                if ioe.errno == 413 and n_frame_of_page > n_frame_in_chunk:
                    n_frame_of_page = max(n_frame_in_chunk, n_frame_of_page // 2)
                    continue
                else:
                    raise IOError(f"Error retrieving data: {ioe.errno}")

            n_frame_in_data = data.shape[0]
            frame_slice = slice(i_frame_in_result, i_frame_in_result + n_frame_in_data)
            result[frame_slice, slice(None), slice(None)] = data

            i_frame_in_result += n_frame_in_data
            i_frame_of_page += n_frame_of_page

            # make sure page starts on a stride step
            offset = i_frame_of_page % stride
            if offset != 0:
                i_frame_of_page += stride - offset

            if i_frame_of_page >= frame_slice.stop:
                break

        return result

    def read(self, n_frames=None, stride=None, atom_indices=None) -> namedtuple:
        """
        Compatibility with mdtraj.formats.hdf5.H5TrajectoryFile.read

        :returns namedtuple:
             The returned namedtuple will have the fields "coordinates", "time", "cell_lengths",
             "cell_angles", "velocities", "kineticEnergy", "potentialEnergy",
             "temperature" and "alchemicalLambda". Each of the fields in the
             returned namedtuple will either be a numpy array or None, dependening
             on if that data was saved in the trajectory. All of the data shall be
             n units of "nanometers", "picoseconds", "kelvin", "degrees" and
             "kilojoules_per_mole".
        """
        if n_frames is None:
            n_frames = numpy.inf
        if stride is not None:
            stride = int(stride)

        total_n_frames = self.get_n_frame()
        frame_slice = slice(0, min(n_frames, total_n_frames), stride)
        if frame_slice.stop - frame_slice.start == 0:
            return []

        if atom_indices is None:
            # get all the atoms
            atom_indices = slice(None)
        # TODO: check atom_indices for type and fit to n_frames

        keys = [field["key"] for field in self.fields]
        Frames = namedtuple("Frames", keys)
        kwargs = {}
        for key in keys:
            if self.has_dataset(key):
                kwargs[key] = self.read_atom_dataset_progressively(
                    key, frame_slice, atom_indices
                )
            else:
                kwargs[key] = None
        return Frames(**kwargs)

    def read_frame_slice_as_traj(
        self, frame_slice, atom_indices, coordinates_only=False
    ):
        """
        :param frame_slice: int | slice(0, n, stride) | [int]
        :param atom_indices: list[int] | None
        """

        kwargs = Dict(topology=self.fetch_topology(atom_indices))
        logger.info(f"{kwargs['topology']}")

        if atom_indices is None:
            atom_indices = slice(None)

        logger.info(tic("reading frames"))
        if coordinates_only:
            fields = [field for field in self.fields if field["key"] == "coordinates"]
        else:
            fields = self.fields

        for field in fields:
            if "traj_key" in field and self.has_dataset(field.key):
                dataset = self.get_dataset(field.key)

                if len(field.shape) == 3:
                    if isinstance(frame_slice, int) or check_is_int_list(frame_slice):
                        data = dataset[frame_slice, atom_indices, :]
                    else:
                        data = self.read_atom_dataset_progressively(
                            field.key, frame_slice, atom_indices
                        )
                elif len(field.shape) == 2:
                    data = dataset[frame_slice, :]
                else:
                    data = dataset[frame_slice]

                kwargs[field.traj_key] = data
        logger.info(toc())

        traj = Trajectory(**kwargs)

        return traj

    def read_as_traj(self, stride=1, atom_indices=None):
        if atom_indices is None:
            atom_indices = self.atom_indices
        n = self.get_n_frame()
        # NOTE: slice can only start on 0
        return self.read_frame_slice_as_traj(slice(0, n, stride), atom_indices)

    def read_frame_as_traj(self, i_frame, atom_indices=None):
        if atom_indices is None:
            atom_indices = self.atom_indices

        if i_frame < 0:
            i_frame = self.get_n_frame() + i_frame

        if self.last_i_frame == i_frame:
            if numpy.array_equal(self.last_atom_indices, atom_indices):
                logger.info(f"same as last frame {i_frame}")
                return self.last_frame

        frame = self.read_frame_slice_as_traj(i_frame, atom_indices)

        self.last_frame = frame
        self.last_i_frame = i_frame
        self.last_atom_indices = atom_indices

        return frame

    def read_parmed_box_vectors(self, i_frame):
        """
        :return: iterable of 3, each 3-element tuple of unit cell vectors in angstroms

        """
        if not self.has_dataset("cell_lengths") or not self.has_dataset("cell_angles"):
            return None
        these_cell_lengths = self.get_dataset("cell_lengths")[i_frame]
        these_cell_angles = self.get_dataset("cell_angles")[i_frame]
        v1, v2, v3 = lengths_and_angles_to_box_vectors(
            these_cell_lengths[0],  # a
            these_cell_lengths[1],  # b
            these_cell_lengths[2],  # c
            these_cell_angles[0],  # alpha
            these_cell_angles[1],  # beta
            these_cell_angles[2],  # gamma
        )
        unitcell_vectors = numpy.swapaxes(numpy.dstack((v1, v2, v3)), 1, 2) * 10
        return unitcell_vectors[0]

    def write_header(self, n_atoms, keys):
        for k, v in self.default_attrs.items():
            self.set_attr(k, v)

        for field in self.fields:
            if field.key not in keys:
                continue

            frame_shape = tuple(field.shape[1:])
            if len(frame_shape) > 1 and frame_shape[0] is None:
                frame_shape = (n_atoms, *frame_shape[1:])

            self.create_extendable_dataset(field.key, frame_shape, numpy.float32)
            self.set_attr("units", field.units, dataset_key=field.key)

    def write(
        self,
        coordinates,
        time=None,
        cell_lengths=None,
        cell_angles=None,
        velocities=None,
        kineticEnergy=None,
        potentialEnergy=None,
        temperature=None,
        alchemicalLambda=None,
    ):
        """
        Compatibility with mdtraj.formats.hdf5.H5TrajectoryFile.write
        """
        logger.info(tic("writing coordinates"))
        frames_by_key = {}
        for field in self.fields:
            # trick to read the args by their string name
            frames = locals().get(field.key)
            if frames is not None:
                frames_by_key[field.key] = ensure_type(
                    in_units_of(frames, None, field.units),
                    name=field.key,
                    dtype=numpy.float32,
                    shape=field.shape,
                    ndim=len(field.shape),
                    add_newaxis_on_deficient_ndim=True,
                    warn_on_cast=False,
                )

        if not self.has_header:
            n_atoms = frames_by_key["coordinates"].shape[-2]
            self.write_header(n_atoms, list(frames_by_key.keys()))
            self.has_header = True

        for key, frames in frames_by_key.items():
            self.extend_dataset(key, frames)

        self.flush()
        logger.info(toc())
        logger.info(f'final shape: {frames_by_key["coordinates"].shape}')


class EasyTrajH5Reporter(_BaseReporter):
    @property
    def backend(self):
        return EasyTrajH5File


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


def get_mdtraj_topology_from_dict(topology_dict) -> Topology:
    topology = Topology()

    for chain_dict in sorted(topology_dict["chains"], key=operator.itemgetter("index")):
        chain = topology.add_chain()
        for residue_dict in sorted(
            chain_dict["residues"], key=operator.itemgetter("index")
        ):
            try:
                ref_seq = residue_dict["resSeq"]
            except KeyError:
                ref_seq = None
                logger.warning(
                    "No resSeq information found in HDF h5, defaulting to zero-based indices"
                )
            try:
                segment_id = residue_dict["segmentID"]
            except KeyError:
                segment_id = ""
            residue = topology.add_residue(
                residue_dict["name"], chain, resSeq=ref_seq, segment_id=segment_id
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


def convert_h5_to_dcd_and_pdb(h5_fname, is_solvent=True):
    dcd = Path(h5_fname).with_suffix(".dcd")
    pdb = Path(h5_fname).with_suffix(".pdb")
    yaml = Path(h5_fname).with_suffix(".rselect.yaml")

    is_convert = False
    meta = dict(is_solvent=is_solvent)
    if yaml.exists():
        old_meta = load_yaml(yaml)
        if meta != old_meta:
            is_convert = True
    if not dcd.exists():
        is_convert = True

    if is_convert:
        logger.info(f"Loading {h5_fname}")

        mask = "not {solvent}" if not is_solvent else ""
        traj = EasyTrajH5File(h5_fname, atom_mask=mask).read_as_traj()

        try:
            logger.info("Imaging frames")
            traj.image_molecules(inplace=True)
        except Exception:
            logger.error("Failed to image molecule")

        try:
            logger.info("Aligning frames")
            traj = traj.superpose(traj[0], atom_indices=traj.top.select("protein"))
        except Exception:
            logger.error("Failed to align")

        logger.info("Centering")
        traj.center_coordinates()

        logger.info(f"Saving {dcd} {pdb}")
        traj.save_dcd(dcd)

        traj[0].save_pdb(pdb)
        remove_model_lines(pdb)
        dump_yaml(meta, yaml)

    return dcd, pdb
