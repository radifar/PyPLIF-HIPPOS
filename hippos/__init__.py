"""
PyPLIF HIPPOS
HIPPOS Is PyPLIF On Steroids. A Molecular Interaction Fingerprinting Tool for Docking Results of Autodock Vina and PLANTS.
"""

# Solution to use simple import for dual way installation
# https://stackoverflow.com/a/49375740/11445093
import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
