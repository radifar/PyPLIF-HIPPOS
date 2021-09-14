"""
Test for functions in pyplif_hippos.initialize.parse_docking_conf
"""

# Import package, test suite, and other packages as needed
import os
import sys

import pytest

try:
    # Open Babel >= 3.0
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob

from pyplif_hippos import get_direct_bitstring, get_bitstring, Residue


def test_get_direct_bitstring(hippos_config):
    """
    Simple test for get_direct_bitstring
    specific test for various molecular format will be carried out in another part.
    """

    # Arrange

    ligand_mol_list = []
    mol_path = "tests/data/direct_ifp/mol2_vina/"
    protein_name = mol_path + "protein_vina.mol2"
    ligand_names = [
        mol_path + "vina1.mol2",
        mol_path + "vina2.mol2",
        mol_path + "vina3.mol2",
        mol_path + "vina4.mol2",
        mol_path + "vina5.mol2",
    ]

    convert = ob.OBConversion()
    convert.SetInFormat("mol2")
    protein_mol = ob.OBMol()
    convert.ReadFile(protein_mol, protein_name)

    for ligand in ligand_names:
        mol = ob.OBMol()
        convert.ReadFile(mol, ligand)
        ligand_mol_list.append(mol)

    # Act

    residues = get_direct_bitstring(protein_mol, ligand_mol_list, hippos_config)

    # Assert

    assert residues["ARG116"].AA_name == "ARG"


def test_residue_class(hippos_config):
    """
    Simple test on Residue class
    """

    # Arrange

    custom_settings = {
        "omit_interaction": hippos_config.omit_interaction,
        "backbone": hippos_config.use_backbone,
        "res_weight": hippos_config.res_weight,
        "output_mode": hippos_config.output_mode,
    }

    mol_path = "tests/data/direct_ifp/mol2_vina/"
    protein_name = mol_path + "protein_vina.mol2"
    convert = ob.OBConversion()
    convert.SetInFormat("mol2")
    protein_mol = ob.OBMol()
    convert.ReadFile(protein_mol, protein_name)

    # Act

    residues = {}
    for name, num, in zip(hippos_config.residue_name, hippos_config.residue_number):
        residues[name] = Residue(protein_mol, name, num, custom_settings)

    # Assert

    assert residues["ARG116"].AA_name == "ARG"
