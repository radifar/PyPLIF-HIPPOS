"""
The tests for omit interaction feature
"""

import os
import sys

from collections import namedtuple

from pyplif_hippos import ParseConfig, hippos, similarity


def test_configuration_single_omit_interaction(tmpdir):
    """Test configuration for omitting specific interaction"""

    # Arrange

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write(
        """
docking_method    plants    # plants or vina
docking_conf      plants-003.conf

similarity_coef   tanimoto mcconnaughey

full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

omit_interaction  hydrophobic  ARG223

full_outfile plants_full_ifp.csv
sim_outfile plants_similarity.csv
logfile plants.log
"""
    )

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    # Act

    hippos_config = ParseConfig()
    hippos_config.parse_config()
    omit_interaction = hippos_config.omit_interaction[0]

    # Assert

    assert omit_interaction.interaction_type == "hydrophobic"
    assert omit_interaction.res_name == ["ARG223"]


def test_configuration_omit_multiple_residue_interaction(tmpdir):
    """Test configuration for omitting multiple residue interaction"""

    # Arrange

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write(
        """
docking_method    plants    # plants or vina
docking_conf      plants-003.conf

similarity_coef   tanimoto mcconnaughey

full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

omit_interaction  hydrophobic  ARG150 TRP177 ARG223

full_outfile plants_full_ifp.csv
sim_outfile plants_similarity.csv
logfile plants.log
"""
    )

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    # Act

    hippos_config = ParseConfig()
    hippos_config.parse_config()
    omit_interaction = hippos_config.omit_interaction[0]

    # Assert

    assert omit_interaction.interaction_type == "hydrophobic"
    assert omit_interaction.res_name == ["ARG150", "TRP177", "ARG223"]


def test_configuration_omit_multiple_interaction_type(tmpdir):
    """Test configuration for omitting multiple interaction type"""

    # Arrange

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write(
        """
docking_method    plants    # plants or vina
docking_conf      plants-003.conf

similarity_coef   tanimoto mcconnaughey

full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

omit_interaction  hydrophobic  ARG223
omit_interaction  h_bond  ARG292

full_outfile plants_full_ifp.csv
sim_outfile plants_similarity.csv
logfile plants.log
"""
    )

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    # Act

    hippos_config = ParseConfig()
    hippos_config.parse_config()
    omit_interaction_1 = hippos_config.omit_interaction[0]
    omit_interaction_2 = hippos_config.omit_interaction[1]

    # Assert

    assert omit_interaction_1.interaction_type == "hydrophobic"
    assert omit_interaction_1.res_name == ["ARG223"]

    assert omit_interaction_2.interaction_type == "h_bond"
    assert omit_interaction_2.res_name == ["ARG292"]


def test_configuration_long_interaction_type(tmpdir):
    """Test configuration checking all long interaction_type"""

    # Arrange

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write(
        """
docking_method    plants    # plants or vina
docking_conf      plants-003.conf

similarity_coef   tanimoto mcconnaughey

full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

omit_interaction  hydrophobic  ARG116
omit_interaction  aromatic  GLU117
omit_interaction  h_bond   LEU132
omit_interaction  electrostatic  LYS148
omit_interaction  h_bond_donor  ASP149
omit_interaction  h_bond_acceptor  ARG150
omit_interaction  electrostatic_positive  ARG154
omit_interaction  electrostatic_negative  TRP177
omit_interaction  aromatic_facetoface  SER178
omit_interaction  aromatic_edgetoface  ILE221

full_outfile plants_full_ifp.csv
sim_outfile plants_similarity.csv
logfile plants.log
"""
    )

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    # Act

    hippos_config = ParseConfig()
    hippos_config.parse_config()
    omit_interaction_1 = hippos_config.omit_interaction[0]
    omit_interaction_2 = hippos_config.omit_interaction[1]
    omit_interaction_3 = hippos_config.omit_interaction[2]
    omit_interaction_4 = hippos_config.omit_interaction[3]
    omit_interaction_5 = hippos_config.omit_interaction[4]
    omit_interaction_6 = hippos_config.omit_interaction[5]
    omit_interaction_7 = hippos_config.omit_interaction[6]
    omit_interaction_8 = hippos_config.omit_interaction[7]
    omit_interaction_9 = hippos_config.omit_interaction[8]
    omit_interaction_10 = hippos_config.omit_interaction[9]

    # Assert

    assert omit_interaction_1.interaction_type == "hydrophobic"
    assert omit_interaction_1.res_name == ["ARG116"]
    assert omit_interaction_2.interaction_type == "aromatic"
    assert omit_interaction_2.res_name == ["GLU117"]
    assert omit_interaction_3.interaction_type == "h_bond"
    assert omit_interaction_3.res_name == ["LEU132"]
    assert omit_interaction_4.interaction_type == "electrostatic"
    assert omit_interaction_4.res_name == ["LYS148"]
    assert omit_interaction_5.interaction_type == "h_bond_donor"
    assert omit_interaction_5.res_name == ["ASP149"]
    assert omit_interaction_6.interaction_type == "h_bond_acceptor"
    assert omit_interaction_6.res_name == ["ARG150"]
    assert omit_interaction_7.interaction_type == "electrostatic_positive"
    assert omit_interaction_7.res_name == ["ARG154"]
    assert omit_interaction_8.interaction_type == "electrostatic_negative"
    assert omit_interaction_8.res_name == ["TRP177"]
    assert omit_interaction_9.interaction_type == "aromatic_facetoface"
    assert omit_interaction_9.res_name == ["SER178"]
    assert omit_interaction_10.interaction_type == "aromatic_edgetoface"
    assert omit_interaction_10.res_name == ["ILE221"]


def test_configuration_short_interaction_type(tmpdir):
    """Test configuration checking all short interaction_type"""

    # Arrange

    config_file = tmpdir.mkdir("sub").join("config.txt")
    config_file.write(
        """
docking_method    plants    # plants or vina
docking_conf      plants-003.conf

similarity_coef   tanimoto mcconnaughey

full_ref  00000100000000000000000000000000000100000000000001000000000000010000001000000000000000000001000000000000000000000000000000101000000000000000000101000000000010000 00010101000000000000000000000000000100000000000001010000000000010000001000000000000010000000000000000000000001011000001000001000000000000000000101000000000000000 00010101000000100000000000000000000100000000000001010100100000010000001000000000000010000001000000000000010000000000100000101010000000000000000001000000000000000

residue_name ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409
residue_number 40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333

omit_interaction  HPB  ARG116
omit_interaction  ARM  GLU117
omit_interaction  HBD   LEU132
omit_interaction  ELE  LYS148
omit_interaction  HBD_DON  ASP149
omit_interaction  HBD_ACC  ARG150
omit_interaction  ELE_POS  ARG154
omit_interaction  ELE_NEG  TRP177
omit_interaction  ARM_F2F  SER178
omit_interaction  ARM_E2F  ILE221

full_outfile plants_full_ifp.csv
sim_outfile plants_similarity.csv
logfile plants.log
"""
    )

    arg = os.path.join(config_file.dirname, config_file.basename)

    if len(sys.argv) > 1:
        sys.argv[1] = arg
    else:
        sys.argv.append(arg)

    # Act

    hippos_config = ParseConfig()
    hippos_config.parse_config()
    omit_interaction_1 = hippos_config.omit_interaction[0]
    omit_interaction_2 = hippos_config.omit_interaction[1]
    omit_interaction_3 = hippos_config.omit_interaction[2]
    omit_interaction_4 = hippos_config.omit_interaction[3]
    omit_interaction_5 = hippos_config.omit_interaction[4]
    omit_interaction_6 = hippos_config.omit_interaction[5]
    omit_interaction_7 = hippos_config.omit_interaction[6]
    omit_interaction_8 = hippos_config.omit_interaction[7]
    omit_interaction_9 = hippos_config.omit_interaction[8]
    omit_interaction_10 = hippos_config.omit_interaction[9]

    # Assert

    assert omit_interaction_1.interaction_type == "hydrophobic"
    assert omit_interaction_1.res_name == ["ARG116"]
    assert omit_interaction_2.interaction_type == "aromatic"
    assert omit_interaction_2.res_name == ["GLU117"]
    assert omit_interaction_3.interaction_type == "h_bond"
    assert omit_interaction_3.res_name == ["LEU132"]
    assert omit_interaction_4.interaction_type == "electrostatic"
    assert omit_interaction_4.res_name == ["LYS148"]
    assert omit_interaction_5.interaction_type == "h_bond_donor"
    assert omit_interaction_5.res_name == ["ASP149"]
    assert omit_interaction_6.interaction_type == "h_bond_acceptor"
    assert omit_interaction_6.res_name == ["ARG150"]
    assert omit_interaction_7.interaction_type == "electrostatic_positive"
    assert omit_interaction_7.res_name == ["ARG154"]
    assert omit_interaction_8.interaction_type == "electrostatic_negative"
    assert omit_interaction_8.res_name == ["TRP177"]
    assert omit_interaction_9.interaction_type == "aromatic_facetoface"
    assert omit_interaction_9.res_name == ["SER178"]
    assert omit_interaction_10.interaction_type == "aromatic_edgetoface"
    assert omit_interaction_10.res_name == ["ILE221"]


def test_replace_bit_char():
    """Test bit replacement function for omitted residue"""

    # Arrange

    bitstring = "1000001"

    omit_hydrophobic = [1, 0, 0, 0, 0, 0, 0]
    omit_aromatic = [0, 1, 1, 0, 0, 0, 0]
    omit_h_bond = [0, 0, 0, 1, 1, 0, 0]
    omit_electrostatic = [0, 0, 0, 0, 0, 1, 1]
    omit_h_bond_donor = [0, 0, 0, 1, 0, 0, 0]
    omit_h_bond_acceptor = [0, 0, 0, 0, 1, 0, 0]
    omit_electrostatic_positive = [0, 0, 0, 0, 0, 1, 0]
    omit_electrostatic_negative = [0, 0, 0, 0, 0, 0, 1]
    omit_aromatic_facetoface = [0, 1, 0, 0, 0, 0, 0]
    omit_aromatic_edgetoface = [0, 0, 1, 0, 0, 0, 0]

    # Act

    bitstring_1 = hippos.replace_bit_char(bitstring, omit_hydrophobic)
    bitstring_2 = hippos.replace_bit_char(bitstring, omit_aromatic)
    bitstring_3 = hippos.replace_bit_char(bitstring, omit_h_bond)
    bitstring_4 = hippos.replace_bit_char(bitstring, omit_electrostatic)
    bitstring_5 = hippos.replace_bit_char(bitstring, omit_h_bond_donor)
    bitstring_6 = hippos.replace_bit_char(bitstring, omit_h_bond_acceptor)
    bitstring_7 = hippos.replace_bit_char(bitstring, omit_electrostatic_positive)
    bitstring_8 = hippos.replace_bit_char(bitstring, omit_electrostatic_negative)
    bitstring_9 = hippos.replace_bit_char(bitstring, omit_aromatic_facetoface)
    bitstring_10 = hippos.replace_bit_char(bitstring, omit_aromatic_edgetoface)

    # Assert

    assert bitstring_1 == "n000001"
    assert bitstring_2 == "1nn0001"
    assert bitstring_3 == "100nn01"
    assert bitstring_4 == "10000nn"
    assert bitstring_5 == "100n001"
    assert bitstring_6 == "1000n01"
    assert bitstring_7 == "10000n1"
    assert bitstring_8 == "100000n"
    assert bitstring_9 == "1n00001"
    assert bitstring_10 == "10n0001"


def test_cleanup_omitted_interaction():
    """Test for bitstring preparation prior to similarity calculation"""

    # Arrange

    refbit = "000001000101"
    tgtbit = "11n00n000011"

    # Act

    clean_refbit, clean_tgtbit = similarity.clean_omitted_interactions(refbit, tgtbit)

    # Assert

    assert clean_refbit == "0000000101"
    assert clean_tgtbit == "1100000011"
