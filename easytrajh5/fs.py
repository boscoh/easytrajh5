import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import h5py as h5
import numpy as np
from addict import Dict
from deepdiff import DeepDiff
from parmed import unit
from pydash import py_
from rich.pretty import pretty_repr
from ruyaml import YAML

logger = logging.getLogger(__name__)

yaml = YAML(typ="safe")
yaml.default_flow_style = False


class Encoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, unit.Quantity):
            value = o._value
            if isinstance(value, np.ndarray):
                value = value.tolist()
            elif isinstance(value, list):
                value = np.array(value).tolist()
            result = {"type": "quantity", "value": value, "unit": str(o.unit)}
            if not hasattr(unit, result["unit"]):
                result["unit_repr"] = repr(o.unit)
            return result
        if isinstance(o, Path):
            return str(o)
        return json.JSONEncoder.default(self, o)


def decode_object(o):
    if isinstance(o, dict):
        result = {}
        for key, item in o.items():
            result[key] = decode_object(item)
        return result
    elif isinstance(o, list):
        return [decode_object(item) for item in o]
    return o


def load_json(f, is_addict=False):
    with open(f) as handle:
        result = decode_object(json.load(handle))
    if is_addict:
        result = Dict(result)
    return result


def dump_json(o, f):
    check_file_dir(f)
    with open(f, "w") as handle:
        json.dump(o, handle, cls=Encoder, default=str)


def load_yaml(f, is_addict=False):
    with open(f) as handle:
        result = decode_object(yaml.load(handle))
    if is_addict:
        result = Dict(result)
    return result


def load_yaml_dict(f):
    return load_yaml(f, is_addict=True)


def dump_yaml(o, f, mode="w"):
    check_file_dir(f)
    o = json.loads(json.dumps(o, cls=Encoder))
    with open(f, mode) as handle:
        yaml.dump(o, handle)


def create_new_dataset_h5(f, value, key):
    if isinstance(value, list):
        shape = [len(value)]
    elif isinstance(value, np.ndarray):
        shape = value.shape
    else:
        shape = []
    f.create_dataset(key, maxshape=(None, *shape), data=np.array([value]))


def dump_meta_h5(file, value, key):
    path = Path(file)
    mode = "a" if path.is_file() else "w"  # if the file exist
    with h5.File(path, mode) as f:
        f.attrs[key] = value


def dump_h5(file, value, key):
    path = Path(file)

    if path.is_file():  # if the file exist
        with h5.File(path, "a") as f:
            if key in f.keys():
                old_size = f[key].shape[0]
                f[key].resize(old_size + 1, axis=0)
                f[key][old_size:] = value
            else:
                create_new_dataset_h5(f, value, key)
    else:
        with h5.File(path, "w") as f:
            create_new_dataset_h5(f, value, key)


def get_yaml_str(o, indent=0):
    o = json.loads(json.dumps(o, cls=Encoder, default=str))
    import textwrap
    from io import StringIO

    handle = StringIO()
    yaml.dump(o, handle)
    s = handle.getvalue()
    if indent:
        s = textwrap.indent(s, " " * indent)
    if s.endswith("\n"):
        s = s[:-1]
    return s


def print_yaml(o, indent=2):
    print(get_yaml_str(o, indent=indent))


def ensure_dir(d):
    Path(d).mkdir(exist_ok=True, parents=True)


def check_file_dir(f):
    Path(f).parent.mkdir(exist_ok=True, parents=True)


def clear_dir(d):
    if Path(d).exists():
        shutil.rmtree(d)
    Path(d).mkdir()


def get_empty_path_str(fname):
    fname = Path(fname)
    if fname.exists():
        fname.unlink()
    return str(fname)


def copy_to_dir(src, dest_dir):
    Path(dest_dir).mkdir(exist_ok=True, parents=True)
    dst = Path(dest_dir) / Path(src).name
    if dst.resolve() != Path(src).resolve():
        shutil.copy(src, dst)


def copy_file(source, target, follow_symlinks=True) -> None:
    # `target` can be a symlink when `source` is a symlink and `source` is already in the required directory
    # (which we can happen e.g., foam)
    source, target = Path(source), Path(target)

    if target.exists():
        # Target can be a symlink or a real file - either way we're fine and don't need to do anything.
        # (assuming we don't want to overwrite)
        return

    # `target` can either:
    # - not exist at all (which is fine)
    # - be a bad symlink because it links to a non-existent file (which is bad)
    # if it's a bad symlink then target.parent.resolve() will get the parent of the bad symlink (which may exist).
    # This can arise through a mistake from trying to create a symlink.
    # So we should just assert `target` isn't a symlink
    if target.is_symlink():
        raise ValueError(f"{target} is a bad symlink (resolves to {target.resolve()}")

    # Creating the directory
    target_dir = target.parent.resolve()
    target_dir.mkdir(exist_ok=True, parents=True)

    # If source is a symlink then new file will not be a link unless follow_symlinks = False (not default)
    shutil.copyfile(source, target, follow_symlinks=follow_symlinks)


def run(cmd, env=None):
    print(">>>>", cmd)
    subprocess.run(cmd, shell=True, check=True, env=env)


def run_with_output(cmd):
    bytes_output = subprocess.check_output(cmd.split())
    text = bytes_output.decode("utf-8", errors="ignore")
    return text


def move_up_sub_dir(top_dir, glob_str):
    for f in Path(top_dir).glob(glob_str):
        this_dir = f.parent
        if not this_dir.exists():
            continue
        for sibling_f in this_dir.glob("*"):
            shutil.copy(sibling_f, this_dir.parent)
        shutil.rmtree(this_dir)


def append_to_key_list(a_dict, key, val):
    if key not in a_dict:
        a_dict[key] = []
    a_dict[key].append(val)


def transfer_glob_str(from_dir, glob_str, to_dir):
    """
    Make copies of files in from_dir/glob_str in the directory to_dir

        >>> transfer_glob_str('source', '**/*.html', 'target')
    """
    from_dir = Path(from_dir)
    to_dir = Path(to_dir)
    for f in from_dir.glob(glob_str):
        rel_dir = f.parent.relative_to(from_dir)
        print(f"Copy {f} -> {to_dir / rel_dir}")
        copy_to_dir(f, to_dir / rel_dir)


def get_checked_path(f):
    result = Path(f)
    if not result.exists():
        logger.error(f"Error: couldn't find {f}")
        sys.exit(1)
    return result


saved_dirs = []


def push_dir(new_dir):
    new_dir = Path(new_dir)
    new_dir.mkdir(exist_ok=True, parents=True)
    save_dir = os.getcwd()
    saved_dirs.append(save_dir)
    os.chdir(new_dir)


def pop_dir():
    if len(saved_dirs):
        save_dir = saved_dirs.pop()
        os.chdir(save_dir)


def get_active_branch_name():
    head_dir = Path(__file__).resolve().parent.parent.parent / ".git" / "HEAD"
    if head_dir.exists():
        with head_dir.open("r") as f:
            content = f.read().splitlines()

        for line in content:
            if line[0:4] == "ref:":
                return line.partition("refs/heads/")[2]
    return ""


def overwrite_yaml(options_yaml, options):
    backup_options_yaml = Path(options_yaml).with_suffix(".bak.yaml")
    while backup_options_yaml.exists():
        backup_options_yaml = backup_options_yaml.with_suffix(".bak.yaml")
    logger.info(f"Updated {options_yaml} -> {backup_options_yaml}")
    os.rename(options_yaml, backup_options_yaml)
    dump_yaml(options, options_yaml)


def is_equal_dict(options1, options2):
    diff = DeepDiff(options2, options1)
    if diff:
        print(diff.pretty())
    return not diff


def load_yaml_from_yaml_url(v):
    tokens = v.replace("yaml://", "").split("/")
    yaml_fname, dict_path = tokens[0:2]
    d = load_yaml_dict(yaml_fname)
    return py_.get(d, dict_path)


def walk_object(o):
    if isinstance(o, dict):
        a_dict = o
        for k, v in a_dict.items():
            if isinstance(v, dict) or isinstance(v, list):
                walk_object(v)
            elif isinstance(v, str):
                if v.startswith("yaml://"):
                    a_dict[k] = load_yaml_from_yaml_url(v)
    elif isinstance(o, list):
        a_list = o
        for i in range(len(a_list)):
            v = a_list[i]
            if isinstance(v, dict) or isinstance(v, list):
                walk_object(v)
            elif isinstance(v, str):
                if v.startswith("yaml://"):
                    a_list[i] = load_yaml_from_yaml_url(v)


def parse_dict_with_yaml_url(o):
    parsed_o = Dict(o)
    walk_object(parsed_o)
    return parsed_o


time_store = []


def tic(msg="No message"):
    ti = time.perf_counter_ns()  # final time
    time_store.append((msg, ti))
    return f"{msg}:started..."


def toc():
    tf = time.perf_counter_ns()  # final time
    if not len(time_store):
        delta_ms = 0
        msg = "no tic"
    else:
        msg, ti = time_store.pop()
        delta_ns = tf - ti
        delta_ms = round(delta_ns / 1000000)
    return f"{msg}:finished in {delta_ms}ms"


def repr_lines(o, prefix=""):
    lines = pretty_repr(o).split("\n")
    lines[0] = prefix + lines[0]
    return lines


def get_time_str(seconds):
    days = seconds // 86400
    seconds -= days * 86400
    hours = seconds // 3600
    seconds -= hours * 3600
    minutes = seconds // 60
    seconds -= minutes * 60
    if days > 0:
        value = "%d:%d:%02d:%02d" % (
            days,
            hours,
            minutes,
            seconds,
        )
    elif hours > 0:
        value = "%d:%02d:%02d" % (
            hours,
            minutes,
            seconds,
        )
    elif minutes > 0:
        value = "%02d:%02d" % (minutes, seconds)
    elif seconds >= 1:
        value = "0:%02d" % seconds
    elif seconds == 0:
        value = "0:00"
    else:
        value = "0:0%.2f" % seconds
    return value
