"""
Test for dictionary in SIMILARITY_FORMULA
"""

# Import package, test suite, and other packages as needed
from pyplif_hippos import sim_dict
import pytest
import os, sys

def test_similarity_formula():
    """Test if key & value format are correct"""

    assert isinstance(sim_dict, dict)
