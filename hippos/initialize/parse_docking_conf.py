import os, sys

try:
    # Open Babel >= 3.0
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob


def parse_vina_conf(vina_conf):
    original_path = os.getcwd()

    protein_file = ""
    ligand_file = ""
    out = ""

    try:
        configread = open(vina_conf, 'r')
        conf_path = os.path.abspath(vina_conf)
        os.chdir(os.path.dirname(conf_path))
    except FileNotFoundError:
        print("The VINA config file: '%s' can not be found" % (vina_conf))
        sys.exit(1)

    config_lines = [line for line in configread]
    for line in config_lines:
        uncommented = line.split('#')[0]
        line_list = uncommented.split()

        if not line_list:
            continue
        option = line_list[0]

        if option == "receptor":
            protein_file = line_list[2]
        if option == "ligand":
            ligand_file = line_list[2]
            out = ligand_file[:-6] + "_out" + ligand_file[-6:]
        if option == "out":
            out = line_list[2]

    scorelist = []

    try:
        ligand_out_lines = open(out, 'r')
    except FileNotFoundError:
        print("Ligand output file: '%s' can not be found" % (out))
        sys.exit(1)

    for line in ligand_out_lines:
        line = line.split()
        if len(line) > 2:
            if line[2] == "RESULT:":
                scorelist.append(line[3])

    convert = ob.OBConversion()
    convert.SetInFormat("pdbqt")

    protein = ob.OBMol()
    ligands = ob.OBMol()

    convert.ReadFile(protein, protein_file)

    docked_ligands = []
    docked_proteins = []

    not_at_end = convert.ReadFile(ligands, out)
    while not_at_end:
        docked_ligands.append(ligands)
        ligands = ob.OBMol()
        not_at_end = convert.Read(ligands)

    mollist = []
    ligand_name = ligand_file[:-6]
    for name in range(len(docked_ligands)):
        ligand_pose = ligand_name + '_' + str(name + 1)
        mollist.append(ligand_pose)

    docking_results = {
        'protein': protein,
        'docked_ligands': docked_ligands,
        'docked_proteins': docked_proteins,
        'mollist': mollist,
        'scorelist': scorelist,
    }

    os.chdir(original_path)
    return docking_results


def parse_plants_conf(plants_conf):
    original_path = os.getcwd()

    protein_file = ""
    ligand_file = ""
    write_multi_mol2 = True
    plants_output = ""

    try:
        configread = open(plants_conf, 'r')
        conf_path = os.path.abspath(plants_conf)
        os.chdir(os.path.dirname(conf_path))
    except FileNotFoundError:
        print("The PLANTS config file: '%s' can not be found" % (plants_conf))
        sys.exit(1)

    config_lines = [line for line in configread]
    for line in config_lines:
        uncommented = line.split('#')[0]
        line_list = uncommented.split()

        if not line_list:
            continue
        option = line_list[0]

        if option == "output_dir":
            plants_output = line_list[1]
        if option == "protein_file":
            protein_file = line_list[1]
        if option == "ligand_file":
            ligand_file = line_list[1]
        if option == "write_multi_mol2":
            if line_list[1] == "0":
                write_multi_mol2 = False

    convert = ob.OBConversion()
    convert.SetInFormat("mol2")

    protein = ob.OBMol()
    ligands = ob.OBMol()
    flexibles = ob.OBMol()

    convert.ReadFile(protein, protein_file)
    protein.DeleteNonPolarHydrogens()

    # Extracting molecule list & score list
    try:
        os.chdir(plants_output)
        ligand_poses = open('features.csv', 'r')

    except FileNotFoundError:
        print('The protein ligand folder can not be found')
        sys.exit(1)

    # just removing the first line
    ligand_poses.readline()

    mollisttemp = [line for line in ligand_poses]
    mollist = []
    scorelist = []

    for mol in mollisttemp:
        mollist.append(mol.split(',')[0])
        scorelist.append(mol.split(',')[1])

    docked_ligands = []
    docked_proteins = []

    if write_multi_mol2:
        not_at_end = convert.ReadFile(ligands, "docked_ligands.mol2")
        while not_at_end:
            ligands.DeleteNonPolarHydrogens()
            docked_ligands.append(ligands)
            ligands = ob.OBMol()
            not_at_end = convert.Read(ligands)

        not_at_end = convert.ReadFile(flexibles, "docked_proteins.mol2")
        while not_at_end:
            docked_proteins.append(flexibles)
            flexibles = ob.OBMol()
            not_at_end = convert.Read(flexibles)

    else:
        for mol in mollist:
            lig_file = mol + ".mol2"
            flex_file = mol + "_protein.mol2"

            ligands = ob.OBMol()
            flexibles = ob.OBMol()

            convert.ReadFile(ligands, lig_file)
            ligands.DeleteNonPolarHydrogens()
            docked_ligands.append(ligands)

            convert.ReadFile(flexibles, flex_file)
            docked_proteins.append(flexibles)

    docking_results = {
        'protein': protein,
        'docked_ligands': docked_ligands,
        'docked_proteins': docked_proteins,
        'mollist': mollist,
        'scorelist': scorelist,
    }

    os.chdir(original_path)
    return docking_results
