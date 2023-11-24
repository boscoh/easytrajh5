import logging
from path import Path

import appdirs
import pydash as py_
from addict import Dict

from easytrajh5.fs import dump_yaml, load_yaml

logger = logging.getLogger(__name__)

DEFAULT = Dict(
    chimera="/Applications/Chimera.app/Contents/MacOS/chimera",
    pymol="pymol",
    vmd="/Applications/VMD\ 1.9.4a51-arm64-Rev9.app/Contents/vmd/vmd_MACOSXARM64",
)

binary_yaml = appdirs.user_data_dir("rseed.binary.yaml")

if Path(binary_yaml).exists():
    logger.info(f"user configs {binary_yaml}")
    binary = load_yaml(binary_yaml)
else:
    binary = py_.clone_deep(DEFAULT)
    logger.info(f"create binary config file {binary_yaml}")
    dump_yaml(binary, binary_yaml)


def get_binary(key, test_binary=None):
    if test_binary is not None:
        if Path(test_binary).exists():
            return test_binary
    result = binary.get(key)
    print(f"Binaries config `{binary_yaml}`")
    print(f"Entry for {key}: `{result}`")
    return result
