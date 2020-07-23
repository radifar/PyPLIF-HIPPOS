"""
Unit and regression test for the pyplif_hippos package.
"""

# Import package, test suite, and other packages as needed
import hippos
import pytest
import sys

def test_pyplif_hippos_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "hippos" in sys.modules
