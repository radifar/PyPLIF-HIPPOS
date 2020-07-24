"""
PyPLIF HIPPOS
HIPPOS Is PyPLIF On Steroids. A Molecular Interaction Fingerprinting Tool for Docking Results of Autodock Vina and PLANTS.
"""

from __future__ import print_function

from openbabel import OBMol, OBConversion

# Solution to use simple import for dual way installation
# https://stackoverflow.com/a/49375740/11445093
import os, sys
from openbabel import OBMol, OBConversion
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from initialize.parse_conf import parse_config, parse_config_genref
from initialize.parse_docking_conf import parse_vina_conf, parse_plants_conf
from SIMILARITY_FORMULA import sim_dict

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
