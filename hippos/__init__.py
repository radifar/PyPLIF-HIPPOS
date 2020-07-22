"""
PyPLIF HIPPOS
HIPPOS Is PyPLIF On Steroids. A Molecular Interaction Fingerprinting Tool for Docking Results of Autodock Vina and PLANTS.
"""

# Add imports here
from .hippos import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
