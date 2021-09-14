"""
Unit and regression test for the pyplif_hippos package.
"""


import pytest
import sys

import pyplif_hippos


def test_pyplif_hippos_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pyplif_hippos" in sys.modules


def test_collect_ligand(hippos_config):
    """Test collect_ligand in pyplif_hippos module"""

    # Arrange

    # Act

    # Assert
    assert hippos_config.use_backbone == True
