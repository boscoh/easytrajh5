import logging

import mdtraj

from .traj import EasyTrajH5File

logger = logging.getLogger(__name__)


class TrajectoryManager:
    """
    Abstraction to manage multiple trajectories that
    have the same topologies. The exact frame is then referenced
    by [i_frame, i_traj, align_atom_mask].
    """

    def __init__(
            self,
            paths: [str],
            mode: str = "a",
            atom_mask: str = "not {solvent}",
            is_dry_cache: bool = True,
    ) -> None:
        self.paths = paths
        self.mode = mode
        self.atom_mask = atom_mask
        self.is_dry_cache = is_dry_cache

        self.traj_file_by_i = {}
        self.last_i_frame_traj = None
        self.last_frame = None

    def get_n_trajectories(self) -> int:
        return len(self.paths)

    @property
    def backend(self):
        return EasyTrajH5File

    def get_traj_file(self, i_traj) -> EasyTrajH5File:
        if i_traj not in self.traj_file_by_i:
            path = self.paths[i_traj]
            traj_file = self.backend(
                path,
                self.mode,
                atom_mask=self.atom_mask,
                is_dry_cache=self.is_dry_cache,
            )
            if not len(traj_file.get_dataset_keys()):
                raise FileNotFoundError(f"h5: {path} is empty")
            self.traj_file_by_i[i_traj] = traj_file
        return self.traj_file_by_i[i_traj]

    def get_n_frame(self, i_traj) -> int:
        return self.get_traj_file(i_traj).get_n_frame()

    def read_as_frame_traj(self, i_frame_traj) -> mdtraj.Trajectory:
        if self.last_i_frame_traj == i_frame_traj:
            logger.info(f"same as last frame {i_frame_traj}")
            return self.last_frame

        i_frame, i_traj = i_frame_traj[:2]
        if i_frame < 0:
            i_frame = self.get_n_frame(i_traj) + i_frame

        frame = self.get_traj_file(i_traj).read_frame_as_traj(i_frame)

        self.last_i_frame_traj = i_frame_traj
        self.last_frame = frame

        return frame

    def close(self):
        for h5 in self.traj_file_by_i.values():
            h5.close()
