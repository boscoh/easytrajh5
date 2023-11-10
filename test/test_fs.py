import uuid

from easytrajh5 import fs
from path import Path

this_dir = Path(__file__).parent


def test_json():
    fname = this_dir / "test.json"
    d = {"a": "b", "c": {"d": "e"}}
    fs.dump_json(d, fname)
    assert fs.is_equal_dict(d, fs.load_json(fname))


def test_yaml():
    fname = this_dir / "test.yaml"
    x = {"a": "b"}
    fs.dump_yaml(x, fname)
    assert fs.is_equal_dict(x, fs.load_yaml(fname))


def test_make_dirs():
    top_dir = this_dir / "a"
    top_dir.rmtree_p()
    assert not top_dir.exists()

    new_dir = top_dir / "b" / "c"
    fs.ensure_dir(new_dir)
    assert new_dir.exists()


def test_clear_dir():
    source_dir = this_dir / "a"
    fs.ensure_dir(source_dir / "b")
    fnames = make_source_files(source_dir)
    fs.clear_dir(source_dir)
    assert not len(source_dir.listdir())


def make_source_files(source_dir):
    fs.clear_dir(source_dir)
    fnames = [f"{uuid.uuid4()}.txt" for i in range(5)]
    for f in fnames:
        (source_dir / f).write_text("text")
    return fnames


def test_copy_file():
    source_dir = this_dir / "a" / "b" / "c"
    fnames = make_source_files(source_dir)
    target_dir = this_dir / "copy_c"
    target_dir.rmtree_p()
    for f in source_dir.glob("*"):
        fs.copy_file(f, target_dir / f.name)
    for f in fnames:
        assert (target_dir / f).exists()


def test_copy_dir():
    source_dir = this_dir / "a" / "b" / "c"
    fnames = make_source_files(source_dir)
    target_dir = this_dir / "copy_c_2"
    target_dir.rmtree_p()
    for f in source_dir.glob("*"):
        fs.copy_to_dir(f, target_dir)
    for f in fnames:
        assert (target_dir / f).exists()


def test_moveup_sub_dir():
    source_dir = this_dir / "a" / "b" / "c"
    fnames = make_source_files(source_dir)
    nested_dir = source_dir / "dummy" / "d"
    nested_dir.makedirs_p()
    fs.move_up_sub_dir(source_dir, "*")
    assert not source_dir.exists()
    for f in fnames:
        assert (source_dir.parent / f).exists()


def test_transfer_glob_str():
    source_dir = this_dir / "a" / "b" / "c"
    fnames = make_source_files(source_dir)
    target_dir = this_dir / "target"
    target_dir.rmtree_p()
    fs.transfer_glob_str(this_dir / "a", '**/*.txt', target_dir)


if __name__ == "__main__":
    x = {"a": "b"}
    y = {"b": "c"}
    fs.is_equal_dict(x, y)
