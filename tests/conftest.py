import pytest


class Config():
    pass


@pytest.fixture(scope="module")
def hippos_config():
    config = Config()
    config.residue_name = "ARG116 GLU117 LEU132 LYS148 ASP149 ARG150 ARG154 TRP177 SER178 ILE221 ARG223 THR224 GLU226 ALA245 HIS273 GLU275 GLU276 ARG292 ASP294 GLY347 ARG374 TRP408 TYR409".split()
    config.residue_number = "40 41 56 72 73 74 78 101 102 145 147 148 150 169 197 199 200 216 218 271 298 332 333".split()
    config.omit_interaction = []
    config.use_backbone = True
    config.output_mode = dict(full=True, full_nobb=False, simplified=False)
    config.res_weight = []

    return config