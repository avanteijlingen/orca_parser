from setuptools import setup

setup(
    name='orca_parser',
    version='0.1.0',    
    description='A module for parse ORCA output files including hessians (.hess) files',
    url='https://github.com/avanteijlingen/ORCA-Parser',
    author='Alexander van Teijlingen',
    author_email='a.vant@linuxmail.org',
    license='Only for use within ModelMole LTD',
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