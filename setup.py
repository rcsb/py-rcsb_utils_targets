# File: setup.py
# Date: 7-Nov-2020
#
# Updates:
#
#
import re

from setuptools import find_packages
from setuptools import setup

packages = []
thisPackage = "rcsb.utils.targets"

with open("rcsb/utils/targets/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

setup(
    name=thisPackage,
    version=version,
    description="RCSB Python Wrapper Module for Target Utilities",
    long_description="See:  README.md",
    author="John Westbrook",
    author_email="john.westbrook@rcsb.org",
    url="https://github.com/rcsb/py-rcsb_utils_targets",
    #
    license="Apache 2.0",
    classifiers=(
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        # "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ),
    entry_points={"console_scripts": []},
    #
    install_requires=["rcsb.utils.io >= 0.93", "rcsb.utils.config >= 0.35", "rcsb.utils.seq >= 0.48", "chembl-webresource-client >= 0.10.2"],
    packages=find_packages(exclude=["rcsb.mock-data", "rcsb.utils.tests-targets", "rcsb.utils.tests-*", "tests.*"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"]
    },
    #
    test_suite="rcsb.utils.tests-targets",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={"dev": ["check-manifest"], "test": ["coverage"]},
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    # This setting for namespace package support -
    zip_safe=False,
)
