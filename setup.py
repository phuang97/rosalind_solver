import os

from setuptools import setup, find_packages

PACKAGE_ATTR_FILE = "rosalind_solver/__init__.py"


def read(path: str):
    cwd = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(cwd, path), 'r') as fp:
        return fp.read()


def get_attribute(attribute: str):
    for line in read(PACKAGE_ATTR_FILE).splitlines():
        if line.startswith(attribute):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise ValueError(f"cannot find {attribute} in {PACKAGE_ATTR_FILE}")


setup(
    name='rsolver',
    version=get_attribute('__version__'),
    description="phuang97's rosalind solver for practice and documentation",
    author=get_attribute('__author__'),
    author_email=get_attribute('__email__'),
    classifiers=['Programming Language :: Python :: 3.9'],
    packages=['stronghold', 'armory'],
    package_dir={'armory': 'rosalind_solver/armory',
                 'stronghold': 'rosalind_solver/stronghold'},
    include_package_data=True,
    package_data={'': ['data/*']}
)
