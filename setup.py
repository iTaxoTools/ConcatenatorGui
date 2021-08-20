"""The setup module for ConcatenatorQt"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
import pathlib

# Get the long description from the README file
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='concatenator-qt',
    version='0.0.1',
    description='A Qt GUI for Concatenator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Patmanidis Stefanos',
    author_email='stefanpatman91@gmail.com',
    package_dir={'': 'src'},
    packages=find_namespace_packages(
        # exclude=('itaxotools.common*',),
        include=('itaxotools*',),
        where='src',
    ),
    python_requires='>=3.6, <4',
    install_requires=[
        'pyside6>=6.1.1',
        ],
    extras_require={
        'dev': ['pyinstaller==5.0.dev0'],
    },
    entry_points={
        'console_scripts': [
            'concatenator-qt = itaxotools.concatenator.gui:main',
        ],
        'pyinstaller40': [
          'hook-dirs = itaxotools.__pyinstaller:get_hook_dirs',
          'tests = itaxotools.__pyinstaller:get_PyInstaller_tests'
        ]
    },
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    include_package_data=True,
)
