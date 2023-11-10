import uuid

from path import Path

from easytrajh5 import fs

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


def test_copy_files():
    source_dir = this_dir / "a" / "b" / "c"
    fs.clear_dir(source_dir)
    assert not len(source_dir.glob("*"))

    fnames = [f"{uuid.uuid4()}.txt" for i in range(5)]
    for f in fnames:
        (source_dir / f).write_text("text")

    target_dir = this_dir / "copy_c"
    target_dir.rmtree_p()

    for f in source_dir.glob("*"):
        fs.copy_file(f, target_dir / f.name)

    for f in fnames:
        assert (target_dir / f).exists()

    nested_dir = source_dir / "dummy" / "d"
    nested_dir.makedirs_p()

    fs.move_up_sub_dir(source_dir, "*")
    assert not source_dir.exists()
    for f in fnames:
        assert (source_dir.parent / f).exists()


if __name__ == "__main__":
    x = {"a": "b"}
    y = {"b": "c"}
    fs.is_equal_dict(x, y)
