"""
Unit and regression test for the pyplif_hippos package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys


def test_pyplif_hippos_imported():
    """Sample test, will always pass so long as import statement worked"""
    import pyplif_hippos
    assert "pyplif_hippos" in sys.modules
