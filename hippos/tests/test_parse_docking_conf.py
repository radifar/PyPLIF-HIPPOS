"""
Test for functions in hippos.initialize.parse_docking_conf
"""

# Import package, test suite, and other packages as needed
from hippos import parse_vina_conf, parse_plants_conf
import pytest
import os, sys
from openbabel import OBMol

def test_parse_vina_conf():
    """Test parsing AutoDock Vina configuration"""
    print(os.getcwd())
    docking_results = parse_vina_conf('hippos/tests/data/vina/vina-001.conf')

    assert isinstance(docking_results['protein'], OBMol)
    assert isinstance(docking_results['docked_ligands'][0], OBMol)
    assert len(docking_results['docked_ligands']) == 20

    assert docking_results['docked_proteins'] == []
    assert docking_results['scorelist'][0] == '-6.8'
    assert docking_results['mollist'][0] == '247120_1'
