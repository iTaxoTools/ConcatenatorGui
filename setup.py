"""The setup module for ConcatenatorGui"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
import pathlib

# Get the long description from the README file
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='concatenator-gui',
    version='0.2.1',
    description='A Qt GUI for Concatenator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Patmanidis Stefanos',
    author_email='stefanpatman91@gmail.com',
    package_dir={'': 'src'},
    packages=find_namespace_packages(
        include=('itaxotools*',),
        where='src',
    ),
    python_requires='>=3.8.6, <4',
    install_requires=[
        'pyside6>=6.1.1',
        'itaxotools-common>=0.3.3',
        'itaxotools-pygblocks>=0.1.0',
        'concatenator>=0.2.2',
        'mafftpy>=0.2.0',
        'fasttreepy>=0.2.1',
        'sequence_bouncer==1.23.1',
        'biopython',
        'clipkit',
    ],
    extras_require={
        'dev': ['pyinstaller>=4.5.1'],
    },
    entry_points={
        'console_scripts': [
            'concatenator-gui = itaxotools.concatenator_gui:run',
        ]
    },
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',  # noqa
        "Development Status :: 5 - Production/Stable",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3 :: Only',
    ],
    include_package_data=True,
)
