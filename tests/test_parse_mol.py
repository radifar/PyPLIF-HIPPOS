"""
Testing parse_mol modules, consists of 4 functions:

parse_ligands(ligand_list: List[string]) -> List[OBMol]
parse_protein(protein: string) -> OBMol

        Mutate input list
enumerate_ligand_files(
    ligand_pose -> List[None],
    ligand_files -> List[string]
) -> None
enumerate_ligand_file_list(
    ligand_pose -> List[None],
    ligand_file_list -> List[string]
) -> None
"""

from pyplif_hippos import parse_ligands, parse_protein, enumerate_ligand_files, enumerate_ligand_file_list


def test_parse_ligands(vina_ligands_mol2):
    # Arrange

    # Act

    ligand_mol_list = parse_ligands(vina_ligands_mol2)

    # Assert

    assert len(ligand_mol_list) == 5
    assert ligand_mol_list[0].NumAtoms() == 31


def test_parse_protein():
    # Arrange

    mol_path = "tests/data/direct_ifp/mol2_vina/"
    protein_name = mol_path + "protein_vina.mol2"

    # Act

    protein_mol = parse_protein(protein_name)
    arg116 = protein_mol.GetResidue(39)

    # Assert

    assert arg116.GetName() == "ARG116"
    assert protein_mol.NumResidues() == 390


def test_enumerate_ligand_files(vina_ligands_mol2):
    # Arrange

    ligand_pose = ["tests/data/direct_ifp/mol2_vina/vina0.mol2"]

    # Act

    enumerate_ligand_files(ligand_pose, vina_ligands_mol2)

    # Assert

    assert len(ligand_pose) == 6
    assert ligand_pose[5][-10:] == "vina5.mol2"


def test_enumerate_ligand_file_list():
    # Arrange

    ligand_pose = ["tests/data/direct_ifp/mol2_vina/vina0.mol2"]

    ligand_file_list = [
        "tests/data/direct_ifp/ligand_mol2_1.lst",
        "tests/data/direct_ifp/ligand_mol2_2.lst",
    ]

    # Act

    enumerate_ligand_file_list(ligand_pose, ligand_file_list)

    # Assert

    assert len(ligand_pose) == 11
    assert ligand_pose[10][-11:] == "vina10.mol2"
