"""
PyPLIF HIPPOS
HIPPOS Is PyPLIF On Steroids. A Molecular Interaction Fingerprinting Tool for Docking Results of Autodock Vina and PLANTS.
"""

from __future__ import print_function

import os
import sys

try:
    # Open Babel >= 3.0
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob

# Solution to use simple import for dual way installation
# https://stackoverflow.com/a/49375740/11445093
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from initialize.parse_conf import ParseBase, ParseConfig, ParseConfigGenref
from initialize.parse_docking_conf import parse_vina_conf, parse_plants_conf
from parse_mol import (
    parse_ligands,
    parse_protein,
    enumerate_ligand_files,
    enumerate_ligand_file_list,
)
from ifp_processing import (
    get_bitstring,
    get_direct_bitstring,
    get_complex_bitstring,
    Residue,
)
from similarity import count_abcdp, how_similar, replace_bit_char
from observer import setup_dict, do_task
from SIMILARITY_FORMULA import sim_dict

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
