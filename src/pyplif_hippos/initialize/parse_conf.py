from __future__ import print_function

import sys

from collections import namedtuple


class ParseBase:
    def __init__(self):
        self.use_backbone = True
        self.residue_name = []
        self.residue_number = []
        self.res_weight = []
        self.config_lines = []
        self.output_mode = dict(full=False, full_nobb=False, simplified=False)

    def read_config(self, config_file, config_type):
        if len(sys.argv) > 1:
            config_file = sys.argv[1]
        else:
            marker = ""
            if config_type == "genref":
                marker = "-GENREF"
            print(
                "HIPPOS%s config file not defined, using %s instead"
                % (marker, config_file)
            )
            print(
                "To change the config file, use it as argument after HIPPOS%s" % marker
            )
            print("Example:")
            print("\t hippos%s <config file>\n" % marker.lower())

        try:
            with open(config_file, "r") as f:
                self.config_lines = [line.split("#")[0].split() for line in f]
        except IOError:
            print("The config file: '%s' can not be found" % config_file)
            sys.exit(1)


class ParseConfig(ParseBase):
    def __init__(self):
        super().__init__()
        self.config = "config.txt"
        self.type = "hippos"
        self.read_config(self.config, self.type)

        self.direct_ifp = False
        self.protein = ""
        self.ligand_files = []
        self.ligand_file_list = []
        self.complex_list = []

        self.docking_score = True
        self.docking_method = ""
        self.docking_conf = ""

        self.similarity_coef = []
        self.simplified_ref = []
        self.full_ref = []
        self.full_nobb_ref = []
        self.omit_interaction = []
        self.Omit_interaction = namedtuple(
            "omit_interaction", "interaction_type res_name"
        )
        self.abbr_interaction = dict(
            HPB="hydrophobic",
            ARM="aromatic",
            HBD="h_bond",
            ELE="electrostatic",
            HBD_DON="h_bond_donor",
            HBD_ACC="h_bond_acceptor",
            ELE_POS="electrostatic_positive",
            ELE_NEG="electrostatic_negative",
            ARM_F2F="aromatic_facetoface",
            ARM_E2F="aromatic_edgetoface",
        )

        self.simplified_outfile = "simplified_ifp.csv"
        self.full_outfile = "full_ifp.csv"
        self.full_nobb_outfile = "full_nobb_ifp.csv"
        self.sim_outfile = "similarity.csv"
        self.logfile = "hippos.log"

        self.direct_ifp_options = [
            "direct_ifp",
            "protein",
            "ligand_files",
            "ligand_list",
            "complex_list",
        ]

        self.docking_options = [
            "docking_score",
            "docking_method",
            "docking_conf",
        ]

        self.fingerprinting_options = [
            "similarity_coef",
            "simplified_ref",
            "full_ref",
            "full_nobb_ref",
            "use_backbone",
            "residue_name",
            "residue_number",
            "omit_interaction",
            "res_weight",
        ]

        self.output_options = [
            "output_mode",
            "simplified_outfile",
            "full_outfile",
            "full_nobb_outfile",
            "sim_outfile",
            "logfile",
        ]

    def parse_config(self):
        for line in self.config_lines:
            if not line:
                continue

            option = line[0]
            single_value = None
            multiple_value = None
            if len(line) > 1:
                single_value = line[1]
                multiple_value = line[1:]

            if option in self.direct_ifp_options:
                if option == "direct_ifp":
                    self.direct_ifp = True
                elif option == "protein":
                    self.protein = single_value
                elif option == "ligand_files":
                    self.ligand_files = multiple_value
                elif option == "ligand_list":
                    self.ligand_file_list = multiple_value
                else:
                    self.complex_list = multiple_value

            elif option in self.docking_options:
                if option == "docking_score":
                    bool_val = ["yes", "no"]
                    if single_value in bool_val:
                        if single_value == "no":
                            self.docking_score = False
                elif option == "docking_method":
                    method = ["vina", "plants"]
                    if single_value in method:
                        self.docking_method = single_value
                    else:
                        print("docking method '%s' is not recognized" % single_value)
                        sys.exit(1)
                else:
                    self.docking_conf = single_value

            elif option in self.fingerprinting_options:
                if option == "similarity_coef":
                    self.similarity_coef = multiple_value
                elif option == "simplified_ref":
                    self.simplified_ref = multiple_value
                elif option == "full_ref":
                    self.full_ref = multiple_value
                elif option == "full_nobb_ref":
                    self.full_nobb_ref = multiple_value
                elif option == "use_backbone":
                    bool_val = ["yes", "no"]
                    if single_value in bool_val:
                        if single_value == "no":
                            self.use_backbone = False
                elif option == "residue_name":
                    self.residue_name = multiple_value
                elif option == "residue_number":
                    self.residue_number = multiple_value
                elif option == "omit_interaction":
                    interaction_type = single_value
                    omitted_residue = line[2:]

                    if interaction_type in self.abbr_interaction:
                        interaction_type = self.abbr_interaction[interaction_type]

                    omit_this = self.Omit_interaction(interaction_type, omitted_residue)
                    self.omit_interaction.append(omit_this)
                else:
                    self.res_weight = multiple_value

            elif option in self.output_options:
                if option == "output_mode":
                    for value in multiple_value:
                        if value in ("full", "full_nobb", "simplified"):
                            self.output_mode[value] = True
                        else:
                            print("output_mode '%s' is not recognized" % value)
                            sys.exit(1)
                elif option == "simplified_outfile":
                    self.simplified_outfile = single_value
                elif option == "full_outfile":
                    self.full_outfile = single_value
                elif option == "full_nobb_outfile":
                    self.full_nobb_outfile = single_value
                elif option == "sim_outfile":
                    self.sim_outfile = single_value
                else:
                    self.logfile = single_value

            else:
                print("Warning: '%s' option is not recognized" % option)

        if not any(self.output_mode.values()):
            self.output_mode["full"] = True


class ParseConfigGenref(ParseBase):
    def __init__(self):
        super().__init__()
        self.config = "genref-config.txt"
        self.type = "genref"
        self.read_config(self.config, self.type)

        self.proteins = ""
        self.ligands = ""

        self.outfile = "genref-results.txt"
        self.logfile = "hippos-genref.log"

        self.input_options = [
            "proteins",
            "ligands",
        ]

        self.fingerprinting_options = [
            "use_backbone",
            "residue_name",
            "residue_number",
            "res_weight",
        ]

        self.output_options = [
            "output_mode",
            "outfile",
            "logfile",
        ]

    def parse_config(self):
        for line in self.config_lines:
            if not line:
                continue

            option = line[0]
            single_value = None
            multiple_value = None
            if len(line) > 1:
                single_value = line[1]
                multiple_value = line[1:]

            if option in self.input_options:
                if option == "proteins":
                    self.proteins = multiple_value
                else:
                    self.ligands = multiple_value

            elif option in self.fingerprinting_options:
                if option == "use_backbone":
                    value = single_value
                    bool_val = ["yes", "no"]
                    if value in bool_val:
                        if value == "no":
                            self.use_backbone = False
                elif option == "residue_name":
                    self.residue_name = multiple_value
                elif option == "residue_number":
                    self.residue_number = multiple_value
                else:
                    self.res_weight = multiple_value

            elif option in self.output_options:
                if option == "output_mode":
                    for value in multiple_value:
                        if value in ("full", "full_nobb", "simplified"):
                            self.output_mode[value] = True
                        else:
                            print("output_mode '%s' is not recognized" % value)
                            sys.exit(1)

                elif option == "outfile":
                    self.outfile = single_value
                else:
                    self.logfile = single_value

            else:
                print("Warning: '%s' option is not recognized" % option)

        if not any(self.output_mode.values()):
            self.output_mode["full"] = True
