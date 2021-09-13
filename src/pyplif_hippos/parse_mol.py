import glob
import re

try:
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob


def parse_ligands(ligand_list):
    ligand_mol_list = []
    file_format = ligand_list[0].split(".")[-1]

    convert = ob.OBConversion()
    convert.SetInFormat(file_format)
    for ligand in ligand_list:
        mol = ob.OBMol()
        convert.ReadFile(mol, ligand)
        ligand_mol_list.append(mol)

    return ligand_mol_list


def parse_protein(protein):
    file_format = protein.split(".")[-1]
    convert = ob.OBConversion()
    convert.SetInFormat(file_format)
    protein_mol = ob.OBMol()
    convert.ReadFile(protein_mol, protein)

    return protein_mol


def enumerate_ligand_files(ligand_pose, ligand_files):
    temp_ligand_list = []
    for ligand in ligand_files:
        for filename in glob.glob(ligand):
            temp_ligand_list.append(filename)
    temp_ligand_list.sort(key=lambda ligand: int(re.sub('\D', '', ligand)))
    ligand_pose.extend(temp_ligand_list)


def enumerate_ligand_file_list(ligand_pose, ligand_file_list):
    for ligand_file in ligand_file_list:
        temp_ligand_list = []
        with open(ligand_file, 'r') as ligands:
            for ligand in ligands:
                ligand = ligand.strip()
                for filename in glob.glob(ligand):
                    temp_ligand_list.append(filename)
            temp_ligand_list.sort(key=lambda ligand: int(re.sub('\D', '', ligand)))
            ligand_pose.extend(temp_ligand_list)
