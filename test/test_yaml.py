from easytrajh5 import fs
from path import Path
from rich.pretty import pprint

this_dir = Path(__file__).parent

def test_yaml():
    fname = this_dir / "test.yaml"
    x = {"a":"b"}
    fs.dump_yaml(x, fname)
    y = fs.load_yaml(fname)
    assert fs.is_equal_dict(x, y)


if __name__ == "__main__":
    x = {"a":"b"}
    y = {"b": "c"}
    fs.is_equal_dict(x, y)
