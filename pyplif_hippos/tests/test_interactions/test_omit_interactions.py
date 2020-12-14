"""
The tests for omit interaction feature
"""

import os
import sys

from pyplif_hippos import parse_config


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

    hippos_config = parse_config()
    omit_interaction = hippos_config["omit_interaction"][0]

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

    hippos_config = parse_config()
    omit_interaction = hippos_config["omit_interaction"][0]

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

    hippos_config = parse_config()
    omit_interaction_1 = hippos_config["omit_interaction"][0]
    omit_interaction_2 = hippos_config["omit_interaction"][1]

    # Assert

    assert omit_interaction_1.interaction_type == "hydrophobic"
    assert omit_interaction_1.res_name == ["ARG223"]

    assert omit_interaction_2.interaction_type == "h_bond"
    assert omit_interaction_2.res_name == ["ARG292"]
