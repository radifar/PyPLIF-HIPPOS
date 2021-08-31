import os
import sys

try:
    # Open Babel >= 3.0
    from openbabel import openbabel as ob
except ImportError:
    import openbabel as ob


def parse_vina_conf(vina_conf):
    try:
        with open(vina_conf, "r") as f:
            for line in f:
                line = line.split("#")[0].split()
                if not line:
                    continue

                option, _, value = line
                if option == "receptor":
                    protein_file = value
                elif option == "ligand":
                    ligand_file = value
                    out = ligand_file[:-6] + "_out" + ligand_file[-6:]
                elif option == "out":
                    out = value
    except IOError:
        print("The VINA config file: '%s' can not be found" % vina_conf)
        sys.exit(1)

    original_path = os.getcwd()
    conf_path = os.path.abspath(vina_conf)
    os.chdir(os.path.dirname(conf_path))

    try:
        scorelist = []
        mollist = []
        ligand_name = ligand_file[:-6]
        pose = 1
        with open(out, "r") as f:
            for line in f:
                line = line.split()
                if (len(line) > 2) and (line[2] == "RESULT:"):
                    ligand_pose = ligand_name + "_" + str(pose)
                    scorelist.append(line[3])
                    mollist.append(ligand_pose)
                    pose += 1
    except IOError:
        print("Ligand output file: '%s' can not be found" % out)
        sys.exit(1)

    protein = ob.OBMol()
    ligands = ob.OBMol()
    convert = ob.OBConversion()
    convert.SetInFormat("pdbqt")
    convert.ReadFile(protein, protein_file)

    docked_ligands = []
    docked_proteins = []

    not_at_end = convert.ReadFile(ligands, out)
    while not_at_end:
        docked_ligands.append(ligands)
        ligands = ob.OBMol()
        not_at_end = convert.Read(ligands)

    docking_results = {
        "protein": protein,
        "docked_ligands": docked_ligands,
        "docked_proteins": docked_proteins,
        "mollist": mollist,
        "scorelist": scorelist,
    }

    os.chdir(original_path)
    return docking_results


def parse_plants_conf(plants_conf):
    original_path = os.getcwd()
    write_multi_mol2 = True

    try:
        with open(plants_conf, "r") as f:
            for line in f:
                uncommented = line.split("#")[0]
                line_list = uncommented.split()

                if not line_list:
                    continue
                line_list = line_list[:2]
                option, value = line_list

                if option == "output_dir":
                    plants_output = value
                if option == "protein_file":
                    protein_file = value
                if option == "ligand_file":
                    ligand_file = value
                if option == "write_multi_mol2":
                    if value == "0":
                        write_multi_mol2 = False
    except IOError:
        print("The PLANTS config file: '%s' can not be found" % plants_conf)
        sys.exit(1)

    conf_path = os.path.abspath(plants_conf)
    os.chdir(os.path.dirname(conf_path))

    convert = ob.OBConversion()
    convert.SetInFormat("mol2")

    protein = ob.OBMol()
    ligands = ob.OBMol()
    flexibles = ob.OBMol()

    convert.ReadFile(protein, protein_file)

    # Extracting molecule list & score list
    try:
        os.chdir(plants_output)
        mollist = []
        scorelist = []
        with open("features.csv", "r") as f:
            f.readline()  # just removing the first line
            for mol in f:
                mol = mol.split(",")
                mollist.append(mol[0])
                scorelist.append(mol[1])
    except IOError:
        print("The protein ligand folder can not be found")
        sys.exit(1)

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
        "protein": protein,
        "docked_ligands": docked_ligands,
        "docked_proteins": docked_proteins,
        "mollist": mollist,
        "scorelist": scorelist,
    }

    os.chdir(original_path)
    return docking_results
