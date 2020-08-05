"""
Test for functions in pyplif_hippos.initialize.parse_docking_conf
"""

# Import package, test suite, and other packages as needed
from pyplif_hippos import parse_vina_conf, parse_plants_conf
import pytest
import os, sys

try:
    # Open Babel >= 3.0
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob

def test_parse_vina_conf():
    """Test parsing AutoDock Vina configuration"""

    docking_results = parse_vina_conf('pyplif_hippos/tests/data/vina/vina-001.conf')

    assert isinstance(docking_results['protein'], ob.OBMol)
    assert isinstance(docking_results['docked_ligands'][0], ob.OBMol)
    assert len(docking_results['docked_ligands']) == 20
    assert docking_results['docked_proteins'] == []

    assert docking_results['scorelist'][0] == '-6.8'
    assert docking_results['mollist'][0] == '247120_1'

    # More test on protein
    protein = docking_results['protein']
    res1 = protein.GetResidue(0)
    assert res1.GetName() == 'GLU'
    

def test_parse_plants_conf():
    """Test parsing PLANTS configuration"""

    docking_results = parse_plants_conf('pyplif_hippos/tests/data/plants/plants-001.conf')

    assert isinstance(docking_results['protein'], ob.OBMol)
    assert isinstance(docking_results['docked_ligands'][0], ob.OBMol)
    assert isinstance(docking_results['docked_proteins'][0], ob.OBMol)
    assert len(docking_results['docked_ligands']) == 25

    assert len(docking_results['docked_proteins']) == 25
    assert docking_results['scorelist'][0] == '-68.2062'
    assert docking_results['mollist'][0] == '247120 CHEMBL344548_entry_00001_conf_01'

    # More test on protein
    protein = docking_results['protein']
    res1 = protein.GetResidue(0)
    assert res1.GetName() == 'GLU77'
