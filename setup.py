from setuptools import setup
exec(open('orca_parser/version.py').read())

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='orca_parser',
    version=__version__,
    description='A module for parse ORCA output files including hessians (.hess) files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/avanteijlingen/ORCA-Parser',
    author='Alexander van Teijlingen',
    author_email='a.vant@linuxmail.org',
    license='BSD 2-clause',
    packages=['orca_parser'],
    install_requires=['ase',
                      'pandas',
                      'numpy',
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.10',
    ],
)
