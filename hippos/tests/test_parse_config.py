"""
Test for functions in hippos.initialize.parse_conf
"""

# Import package, test suite, and other packages as needed
from hippos import parse_config, parse_config_genref
import pytest
import os, sys


def test_parse_config(tmpdir):
    """Test parsing HIPPOS configuration"""

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write("""
docking_method    plants    # plants or vina
docking_conf      plants-003.conf

similarity_coef   tanimoto mcconnaughey

full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

full_outfile plants_full_ifp.csv
sim_outfile plants_similarity.csv
logfile plants.log
""")

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    hippos_config = parse_config()

    assert hippos_config['docking_method'] == 'plants'
    assert hippos_config['docking_conf'] == 'plants-003.conf'
    assert 'tanimoto' in hippos_config['similarity_coef']
    assert 'mcconnaughey' in hippos_config['similarity_coef']

    full_ref1 = '00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000'
    full_ref2 = '00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000'
    full_ref3 = '00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000'
    assert full_ref1 in hippos_config['full_ref']
    assert full_ref2 in hippos_config['full_ref']
    assert full_ref3 in hippos_config['full_ref']
    assert hippos_config['full_nobb_ref'] == []
    assert hippos_config['simplified_ref'] == []
    assert hippos_config['use_backbone'] == True

    residue_name = 'ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409'
    residue_number = '40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333'
    residue_name = residue_name.split()
    residue_number = residue_number.split()
    assert hippos_config['residue_name'] == residue_name
    assert hippos_config['residue_number'] == residue_number

    assert hippos_config['output_mode']['full'] == True
    assert hippos_config['output_mode']['full_nobb'] == False
    assert hippos_config['output_mode']['simplified'] == False

    assert hippos_config['simplified_outfile'] == 'simplified_ifp.csv'
    assert hippos_config['full_outfile'] == 'plants_full_ifp.csv'
    assert hippos_config['full_nobb_outfile'] == 'full_nobb_ifp.csv'

    assert hippos_config['sim_outfile'] == 'plants_similarity.csv'
    assert hippos_config['logfile'] == 'plants.log'


def test_parse_config_genref(tmpdir):
    """Test parsing HIPPOS-genref configuration"""

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write("""
# first residue is 77
residue_name  ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number  40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

proteins    1b9s/protein.mol2 1b9t/protein.mol2 1b9v/protein.mol2
ligands     1b9s/ligand_FDI468_0.mol2 1b9t/ligand_RAI468_0.mol2 1b9v/ligand_RA2468_0.mol2

outfile     ref-results.txt
""")

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    hippos_config = parse_config_genref()

    assert '1b9s/protein.mol2' in hippos_config['proteins']
    assert '1b9t/protein.mol2' in hippos_config['proteins']
    assert '1b9v/protein.mol2' in hippos_config['proteins']

    assert '1b9s/ligand_FDI468_0.mol2' in hippos_config['ligands']
    assert '1b9t/ligand_RAI468_0.mol2' in hippos_config['ligands']
    assert '1b9v/ligand_RA2468_0.mol2' in hippos_config['ligands']

    assert hippos_config['use_backbone'] == True

    residue_name = 'ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409'
    residue_number = '40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333'
    residue_name = residue_name.split()
    residue_number = residue_number.split()
    assert hippos_config['residue_name'] == residue_name
    assert hippos_config['residue_number'] == residue_number

    assert hippos_config['output_mode']['full'] == True
    assert hippos_config['output_mode']['full_nobb'] == False
    assert hippos_config['output_mode']['simplified'] == False

    assert hippos_config['outfile'] == 'ref-results.txt'
    assert hippos_config['logfile'] == 'hippos-genref.log'