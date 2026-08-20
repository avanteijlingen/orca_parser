from setuptools import setup

exec(open("orca_parser/version.py").read())

from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="orca_parser",
    version=__version__,
    description="A module for parsing QM software (originally ORCA) output files including hessians files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/avanteijlingen/orca_parser",
    author="Alexander van Teijlingen, Ross Urquhart",
    author_email="alexvant-public@protonmail.ch",
    license="BSD 2-clause",
    packages=["orca_parser"],
    install_requires=[
        "ase",
        "pandas",
        "numpy",
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.10",
    ],
)
