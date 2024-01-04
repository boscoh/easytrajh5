import json
import re
from typing import Union, KeysView, Any

import h5py
import numpy
import numpy as np
from addict import Dict
from path import Path
from rich import box
from rich.console import Console
from rich.pretty import pprint
from rich.table import Table

from .select import parse_number_list


def convert_str_to_bytes(s):
    if not isinstance(s, bytes):
        s = s.encode("ascii")
    return s


def convert_bytes_to_str(b: Union[numpy.bytes_, bytes, str]):
    if isinstance(b, numpy.bytes_):
        b = b.tobytes()
    if not isinstance(b, str):
        b = b.decode()
    return b


class EasyH5File:
    """
    EasyH5File is a convenience class to make h5py.File objects
    easier to use for the following use-cases:

        - easily create and add to fixed/extendable numpy arrays
        - add strings as byte column dataset
        - add JSON-compatible object literals
        - add binary blobs; which allows insertion of binary files
        - handle attributes/dataset access using method-based lookups
        - add diagnostics such as JSON schema and print function

    More importantly, this object allows the internal handle type
    of h5py.File to be easily overrideable with h5pyd.File for
    processing with an HSDS data server.

    The h5py.File is stored as self.handle.

    Note: remember to .flush() or .close() object to force write
    """

    def __init__(self, fname, mode="a"):
        self.fname = fname
        self.mode = mode
        self.handle = self.open(self.fname, mode=self.mode)
        self.is_open = True

    def open(self, fname, mode):
        # We use this method to open the h5py object
        # to allow this method to be overriable
        # with a reference to h5pyd for remote
        # access with an HSDS server
        return h5py.File(fname, mode)

    def flush(self):
        self.handle.flush()

    def close(self):
        if self.is_open:
            self.handle.close()
            self.is_open = False

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def __del__(self):
        self.close()

    def has_dataset(self, key) -> bool:
        return key in self.handle

    def get_dataset_keys(self) -> [str]:
        def get_node_keys(node_key=None):
            result = []
            if node_key is None:
                node = self.handle
            else:
                node = self.handle[node_key]
            for key in node.keys():
                if node_key is not None:
                    leaf_key = f"{node_key}/{key}"
                else:
                    leaf_key = key
                if isinstance(self.handle[leaf_key], h5py.Group):
                    result.extend(get_node_keys(leaf_key))
                else:
                    result.append(leaf_key)
            return result

        return get_node_keys()

    def set_array_dataset(self, key, array_data):
        self.clear_dataset(key)
        self.handle.create_dataset(key, data=array_data)

    def clear_dataset(self, key):
        if key in self.handle:
            del self.handle[key]

    def get_dataset(self, key) -> h5py.Dataset:
        return self.handle[key]

    def delete_dataset(self, key):
        if self.has_dataset(key):
            del self.handle[key]

    def create_extendable_dataset(self, key, frame_shape=(), dtype=numpy.float32):
        self.handle.create_dataset(
            key,
            chunks=True,
            compression="gzip",
            dtype=dtype,
            shape=(0, *frame_shape),
            maxshape=(None, *frame_shape),
        )

    def extend_dataset(self, key, new_frames):
        shape = self.handle[key].shape
        n_frame = shape[0]
        frame_shape = shape[1:]
        new_frames = numpy.array(new_frames)
        if new_frames.shape == frame_shape:
            self.handle[key].resize(n_frame + 1, axis=0)
            self.handle[key][-1] = new_frames
        elif new_frames.shape[1:] == frame_shape:
            self.handle[key].resize(n_frame + len(new_frames), axis=0)
            self.handle[key][n_frame:] = new_frames
        else:
            raise ValueError(
                f"Shape of new frames {new_frames.shape} doesn't match {frame_shape}"
            )

    def set_bytes_dataset(self, key, blob):
        self.clear_dataset(key)
        self.handle.create_dataset(key, data=numpy.array([blob]))

    def get_bytes_dataset(self, key) -> bytes:
        return self.get_dataset(key)[0].tobytes()

    def insert_file_to_dataset(self, key, fname):
        with open(fname, "rb") as f:
            blob = f.read()
        self.set_bytes_dataset(key, blob)

    def extract_file_from_dataset(self, key, fname):
        blob = self.get_bytes_dataset(key)
        with open(fname, "wb") as f:
            f.write(blob)

    def set_str_dataset(self, key, s):
        self.set_bytes_dataset(key, convert_str_to_bytes(s))

    def get_str_dataset(self, key) -> str:
        return convert_bytes_to_str(self.get_bytes_dataset(key))

    def set_json_dataset(self, key, obj):
        self.set_str_dataset(key, json.dumps(obj, default=str))

    def get_json_dataset(self, key) -> Any:
        return json.loads(self.get_str_dataset(key))

    def get_attrs(self, dataset_key=None) -> h5py.AttributeManager:
        """Returns h5py attrs object"""
        return self.get_dataset(dataset_key).attrs if dataset_key else self.handle.attrs

    def get_attr_keys(self, dataset_key=None) -> KeysView:
        return self.get_attrs(dataset_key).keys()

    def has_attr(self, key, dataset_key=None) -> bool:
        return key in self.get_attrs(dataset_key)

    def get_attr(self, key, dataset_key=None) -> Any:
        attrs = self.get_attrs(dataset_key)
        v = attrs[key]
        if isinstance(v, numpy.bytes_):
            v = convert_bytes_to_str(v)
        if isinstance(v, h5py._hl.base.Empty):
            v = None
        if isinstance(v, numpy.integer):
            v = int(v)
        if hasattr(numpy, "float") and isinstance(v, numpy.float):
            v = float(v)
        return v

    def get_attr_dict(self, dataset_key=None) -> dict:
        """Returns h5py attrs as a dict"""
        return {
            key: self.get_attr(key, dataset_key)
            for key in self.get_attr_keys(dataset_key)
        }

    def set_attr(self, key, value, dataset_key=None):
        attrs = self.get_attrs(dataset_key)
        if key in attrs:
            attrs.modify(key, value)
        else:
            attrs.create(key, value)

    def get_dataset_schema(self, key) -> Dict:
        dataset = self.get_dataset(key)
        config = Dict(key=key)
        config.shape = list(dataset.shape)
        if dataset.chunks is not None:
            config.chunks = list(dataset.chunks)
        if dataset.maxshape[0] is None:
            config.is_extensible = True
            config.frame_shape = list(dataset.maxshape[1:])
            config.n_frame = config.shape[0]
        config.dtype = str(dataset.dtype)
        if re.search(r"\|S\d+", config.dtype):
            config.dtype = f"string({config.dtype[2:]})"
        for attr_key in self.get_attr_keys(key):
            config.attr[attr_key] = self.get_attr(attr_key, key)
        return config

    def get_schema(self) -> Dict:
        structure = Dict(datasets=[])
        for key in self.get_dataset_keys():
            structure.datasets.append(self.get_dataset_schema(key))
        for key in self.get_attr_keys():
            structure.attr[key] = self.get_attr(key)
        return structure.to_dict()

    def print_schema(self):
        pprint(self.get_schema())

    def print_dataset_table(self, title=None):
        table = Table(title=title, box=box.SIMPLE)
        table.add_column("dataset")
        table.add_column("shape")
        table.add_column("dtype")
        table.add_column("size (MB)", justify="right")

        def to_mb(n):
            mb = n / 1024**2
            if mb < 0.01:
                return "<1 KB"
            return f"{mb:.2f} MB"

        total = 0
        for key in self.get_dataset_keys():
            dataset = self.get_dataset(key)
            n_byte = dataset.nbytes
            total += n_byte
            table.add_row(key, str(dataset.shape), str(dataset.dtype), to_mb(n_byte))

        table.add_row()
        table.add_row("total", "", "", to_mb(total))

        print()
        console = Console()
        console.print(table)

    def print_dataset(self, key, frames=None, is_json=True):
        dataset = self.get_dataset(key)
        dtype = str(dataset.dtype)

        print(f"     dataset={key}")
        print(f"     dtype={dtype}")

        is_string = re.search(r"\|S\d+", dtype)
        if is_string:
            str_len = int(dtype[2:])
            print()
            if is_json:
                print("---- JSON --------------------")
                pprint(self.get_json_dataset(key))
            else:
                b = self.get_bytes_dataset(key)
                try:
                    s = convert_bytes_to_str(b)
                    print("---- STRING --------------------")
                    if str_len > 2000:
                        print(s[:1000])
                        print("\n....\n")
                        print(s[-1000:])
                    else:
                        print(s)
                except Exception:
                    import sys

                    print("---- BYTES --------------------")
                    if str_len > 2000:
                        sys.stdout.buffer.write(b[:1000])
                        print("\n....")
                        sys.stdout.buffer.write(b[-1000:])
                    else:
                        sys.stdout.buffer.write(b)
        else:
            print(f"     shape={dataset.shape}")
            print()
            if frames is not None:
                print(f"frames({frames})=")
                i_frames = parse_number_list(frames)
                chunk = dataset[i_frames]
            else:
                chunk = dataset[:]
            print(chunk)


def dump_attr_to_h5(h5_fname, value, key):
    path = Path(h5_fname)
    mode = "a" if path.isfile() else "w"  # if the h5 exist
    with h5py.File(path, mode) as f:
        f.attrs[key] = value


def create_dataset_in_h5_file_with_value(h5_file, value, key):
    if isinstance(value, list):
        shape = [len(value)]
    elif isinstance(value, np.ndarray):
        shape = value.shape
    else:
        shape = []
    h5_file.create_dataset(key, maxshape=(None, *shape), data=np.array([value]))


def dump_value_to_h5(h5_fname, value, key):
    """
    Convenience function to quickly append values to a dataset to an h5 file.
    The file will be created if it doesn't exist. Similarly, with the dataset.
    The file will be opened and then closed so that there are no dangling
    handles.
    """
    path = Path(h5_fname)

    if path.isfile():  # if the h5 exist
        with h5py.File(path, "a") as f:
            if key in f.keys():
                old_size = f[key].shape[0]
                f[key].resize(old_size + 1, axis=0)
                f[key][old_size:] = value
            else:
                create_dataset_in_h5_file_with_value(f, value, key)
    else:
        with h5py.File(path, "w") as f:
            create_dataset_in_h5_file_with_value(f, value, key)
