import setuptools

setuptools.setup(
    name="easytrajh5",
    version="0.1",
    author="Bosco Ho",
    author_email="bosco@redesignscience.com",
    description="Redesign Science Cloud Job Client in Python",
    install_requires=[
        "addict",
        "deepdiff",
        "h5py",
        "mdtraj",
        "numpy",
        "orjson",
        "ParmEd",
        "pydash",
        "rich",
        "ruyaml"

    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",
    include_package_data=True,
    package_data={
        "easytrajh5": [
            "easytrajh5/data/*",
        ]
    },
)
