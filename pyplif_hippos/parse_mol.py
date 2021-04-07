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
