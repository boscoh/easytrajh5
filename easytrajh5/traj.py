import logging
from typing import TypeVar, Sequence, List, Union, Tuple

import numpy
import numpy as np
from addict import Dict
from mdtraj import Topology, Trajectory
from mdtraj.reporters.basereporter import _BaseReporter
from mdtraj.utils import ensure_type, in_units_of

from .fs import tic, toc
from .h5 import EasyH5
from .struct import (
    slice_mdtraj_topology,
    get_dict_from_mdtraj_topology,
    get_mdtraj_topology_from_dict,
)

logger = logging.getLogger(__name__)

_Slice = TypeVar("_Slice", bound=Sequence)


class EasyTrajH5File(EasyH5):
    """
    Interface to read/write an mdtraj h5 file:
       1. Interface to read from h5 as Trajectory:
            - get_traj
            - get_frame_traj
            - get_traj_with_frames
       2. Interface to write to h5 from a backend of openmm._BaseReporter:
            - __init__(file, mode)
            - @property.setter topology
            - distance_unit
            - flush (optional)
            - close
            - write
    atom_mask - selection language for atom selection
    dry_cache - cache dry atoms for fast access without waters
    """

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
        # reversed the mdtraj choice to rename this `lambda` in the file
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
        logger.info(f"{fname=} {mode=} {atom_mask=} {is_dry_cache=}")
        logger.info(tic("opening"))
        super().__init__(fname, mode)
        logger.info(toc())

        # to reach here, we've succesfully loaded an .h5 file
        if mode == "w":
            self.has_header = False
        elif mode in ["a", "r"]:
            self.has_header = self.has_dataset("coordinates")
        else:
            raise ValueError("mode must be one of ['r', 'w', 'a']")

        self._topology = None
        self.topology_by_atom_hash = {}

        self.atom_mask = atom_mask
        self.atom_indices = None
        self.last_i_frame = None
        self.last_frame = None
        self.last_atom_indices = None

        if is_dry_cache:
            if atom_mask and atom_mask != "not {solvent}":
                raise ValueError("Can't set is_dry_cache=True and atom_mask at the same time")
            atom_mask = "not {solvent}"
            if mode == "r":
                logger.info("Warning: can't use dry_atom_cache in read-only mode")
            else:
                dry_atom_indices = self.load_dry_topology_from_cache()
                self.atom_indices = dry_atom_indices
                return

        if atom_mask:
            # By referencing self.topology, the topology gets loaded
            sliced_top, atom_indices = slice_mdtraj_topology(self.topology, atom_mask)
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
            logger.info(tic(f"loading topology from '{self.fname}'"))
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
            logger.info(tic("reading dry-atoms"))
            dry_atom_indices = self.get_dataset("dry_atoms")[:]
            logger.info(toc())

            if len(dry_atom_indices) == 0:
                logger.info("dry_atoms=[] -> no dry_topology")
                dry_atom_indices = None
            elif self.has_dataset("dry_topology"):
                logger.info(tic("loading cached dry topology"))
                topology_dict = self.get_json_dataset("dry_topology")
                dry_top = get_mdtraj_topology_from_dict(topology_dict)
                logger.info(toc())
                atom_hash = tuple(dry_atom_indices)
                self.topology_by_atom_hash[atom_hash] = dry_top
        else:
            # haven't tried calculating dry_topology yet
            dry_top, dry_atom_indices = slice_mdtraj_topology(
                self.topology, "not {solvent}"
            )
            n_atom = sum(1 for _ in self.topology.atoms)
            if len(dry_atom_indices) == n_atom:
                logger.info(tic("no solvent => save dry_atoms=[]"))
                self.set_array_dataset("dry_atoms", np.array([]))
                logger.info(toc())
                dry_atom_indices = None
            else:
                logger.info(tic("saving cached dry topology"))
                self.set_json_dataset(
                    "dry_topology", get_dict_from_mdtraj_topology(dry_top)
                )
                self.set_array_dataset("dry_atoms", np.array(dry_atom_indices))
                logger.info(toc())
                atom_hash = tuple(dry_atom_indices)
                self.topology_by_atom_hash[atom_hash] = dry_top

        return dry_atom_indices

    @property
    def topology(self) -> Topology:
        if self._topology is None:
            self.fetch_topology()
        return self._topology

    @topology.setter
    def topology(self, topology):
        self._topology = topology
        self.set_json_dataset("topology", get_dict_from_mdtraj_topology(topology))
        self.handle.flush()

    def read_atom_dataset_progressively(
            self, key, frame_slice, atom_indices=slice(None)
    ):
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

    def iterate_chunks(
            self,
            chunk: int,
            atom_indices: List[int] = None,
            stride: int = 1,
            coordinates_only: bool = False,
            return_slice: bool = False,
    ) -> Union[Trajectory, Tuple[Trajectory, _Slice]]:
        """
        A generator over the trajectory which yields a copy of the trajectory with n_frames <= chunk. i.e., will serve up the following slices:
        [slice(0, chunk, stride), slice(chunk, 2*chunk, stride)...].

        chunk: int - maximimum size of the chunk. Note: will be modified by `stride` and if `n_frames` is not an exact multiple of `chunk` in the final chunk.
        atom_indices: List[int] - atom indices to return.
        stride: int - the stride of the chunk. Must be divisor of `chunk`.
        coordinates_only: bool - whether to return only the coordinates.  This will save time but may affect your analysis routines.
        return_slice: bool - whether to return the slice used.
        """
        n_frames = self.get_n_frame()
        chunk = min(n_frames, chunk)

        assert isinstance(chunk, int) and isinstance(
            stride, int
        ), "chunk and stride should be integer"
        assert (chunk >= 1) and (stride >= 1), "chunk and stride should be 1 or more"
        assert chunk > stride, "chunk should be greater than stride"
        assert chunk % stride == 0, "chunk should be multiple of stride"

        n_chunks = n_frames // chunk
        if n_chunks * chunk < n_frames:
            n_chunks += 1

        start = 0

        for _ in range(n_chunks):
            end = min(start + chunk, n_frames)
            frame_slice = slice(start, end, stride)
            if return_slice:
                yield (
                    self.read_frame_slice_as_traj(
                        frame_slice, atom_indices, coordinates_only
                    ),
                    frame_slice,
                )
            else:
                yield self.read_frame_slice_as_traj(
                    frame_slice, atom_indices, coordinates_only
                )
            start = end

    def read_frame_slice_as_traj(self, frame_slice, atom_indices, coordinates_only=False):
        """
        :param frame_slice: int | slice(0, n, stride)
        :param atom_indices: list[int] | None
        """

        kwargs = Dict(topology=self.fetch_topology(atom_indices))

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
                    if isinstance(frame_slice, int):
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

        if self.last_i_frame == i_frame and self.last_atom_indices == atom_indices:
            logger.info(f"same as last frame {i_frame}")
            return self.last_frame

        frame = self.read_frame_slice_as_traj(i_frame, atom_indices)

        self.last_frame = frame
        self.last_i_frame = i_frame
        self.last_atom_indices = atom_indices

        return frame

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
        logger.info(tic("writing coordinates"))
        frames_by_key = {}
        for field in self.fields:
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

    def get_n_frame(self):
        return self.get_dataset("coordinates").shape[0]

    def __len__(self):
        return self.get_n_frame()


class EasyTrajH5Reporter(_BaseReporter):
    @property
    def backend(self):
        return EasyTrajH5File