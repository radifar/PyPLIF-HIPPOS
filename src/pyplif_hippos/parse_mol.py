import glob
import sys
import re
from typing import Type, List, Tuple

try:
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob


def parse_ligands(ligand_list: str) -> List[Type[ob.OBMol]]:
    """Extract the ligand file names and create OBMol object for each file.

    Parameters
    ----------
    ligand_list: str
        a string that can be divided into list of string where each string
        is ligand file name

    Returns
    -------
    list
        a list of OBMol object containing molecule from each ligand file

    """

    ligand_mol_list = []
    file_format = ligand_list[0].split(".")[-1]
    if file_format not in "pdb mol2 pdbqt".split():
        print("Can not read from this format.\nAccepted file formats: pdb, pdbqt, and mol2")
        sys.exit(1)

    convert = ob.OBConversion()
    convert.SetInAndOutFormats(file_format, "pdb")
    if file_format == "pdb":
        for ligand in ligand_list:
            mol = ob.OBMol()
            convert.ReadFile(mol, ligand)
            mol.AddHydrogens(False, True, 7.4)
            ligand_mol_list.append(mol)
    else:
        read_to_pdb = ob.OBConversion()
        read_to_pdb.SetInFormat("pdb")
        for ligand in ligand_list:
            mol = ob.OBMol()
            convert.ReadFile(mol, ligand)
            mol.AddHydrogens(False, True, 7.4)
            pdb_string = convert.WriteString(mol)
            mol_pdb = ob.OBMol()
            read_to_pdb.ReadString(mol_pdb, pdb_string)
            ligand_mol_list.append(mol_pdb)
    

    return ligand_mol_list


def parse_protein(protein: str) -> Type[ob.OBMol]:
    """Take protein file name and return an OBMol object for that protein

    Parameters
    ----------
    protein : str
        a protein file name

    Returns
    -------
    ob.OBMol
        OBMol object for protein
    """

    file_format = protein.split(".")[-1]
    if file_format not in "pdb mol2 pdbqt".split():
        print("Can not read from this format.\nAccepted file formats: pdb, pdbqt, and mol2")
        sys.exit(1)
    
    convert = ob.OBConversion()
    convert.SetInAndOutFormats(file_format, "pdb")
    if file_format == "pdb":
        protein_mol = ob.OBMol()
        convert.ReadFile(protein_mol, protein)
        protein_mol.AddHydrogens(False, True, 7.4)
    else:
        read_to_pdb = ob.OBConversion()
        read_to_pdb.SetInFormat("pdb")

        protein_mol_origin = ob.OBMol()
        convert.ReadFile(protein_mol_origin, protein)
        if file_format == "mol2":
            protein_mol_origin.AddHydrogens(False, True, 7.4)
        protein_mol_string = convert.WriteString(protein_mol_origin)
        protein_mol = ob.OBMol()
        read_to_pdb.ReadString(protein_mol, protein_mol_string)
    return protein_mol


def parse_complex(complex_list: List[str]) -> List[Tuple[ob.OBMol]]:
    complex_mol_list = []
    for complex in complex_list:
        protein = parse_protein(complex.split()[0])
        ligand_name = complex.split()[1]
        ligand = parse_ligands([ligand_name])[0]
        protein_ligand = (protein, ligand)
        complex_mol_list.append(protein_ligand)
    
    return complex_mol_list


def enumerate_ligand_files(ligand_pose: List[str], ligand_files: List[str]) -> None:
    """Take ligand file name in ligand_files and add it to the ligand_pose list

    Parameters
    ----------
    ligand_pose : List[str]
        list of ligand poses that will be extended
    ligand_files : List[str]
        list of ligand file name
    """

    temp_ligand_list = []
    for ligand in ligand_files:
        for filename in glob.glob(ligand):
            temp_ligand_list.append(filename)
    temp_ligand_list.sort(key=lambda ligand: int(re.sub("\D", "", ligand)))
    ligand_pose.extend(temp_ligand_list)


def enumerate_ligand_file_list(
    ligand_pose: List[str], ligand_file_list: List[str]
) -> None:
    """Take list of ligand_file_list, extract ligand file name from each of them and add it to the ligand_pose

    Parameters
    ----------
    ligand_pose : List[str]
        list of ligand poses that will be extended
    ligand_file_list : List[str]
        list of ligand_file_list, each of them contain several ligand file name
    """

    for ligand_file in ligand_file_list:
        temp_ligand_list = []
        with open(ligand_file, "r") as ligands:
            for ligand in ligands:
                ligand = ligand.strip()
                for filename in glob.glob(ligand):
                    temp_ligand_list.append(filename)
            temp_ligand_list.sort(key=lambda ligand: int(re.sub("\D", "", ligand)))
            ligand_pose.extend(temp_ligand_list)
