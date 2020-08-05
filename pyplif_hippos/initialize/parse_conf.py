from __future__ import print_function

import sys


def parse_config():
    #  Default configuration
    config = "config.txt"

    #  Identifiers initialization
    docking_method = ""
    docking_conf = ""
    similarity_coef = []
    simplified_ref = []
    full_ref = []
    full_nobb_ref = []

    use_backbone = True
    res_name = []
    res_num = []
    res_weight1 = []
    res_weight2 = []
    res_weight3 = []
    res_weight4 = []
    res_weight5 = []
    docking_score = True

    output_mode = {"full": False, "full_nobb": False, "simplified": False}
    output_mode_undefined = {"full": False, "full_nobb": False, "simplified": False}
    '''
    Default values
    '''

    simplified_outfile = "simplified_ifp.csv"
    full_outfile = "full_ifp.csv"
    full_nobb_outfile = "full_nobb_ifp.csv"
    sim_outfile = "similarity.csv"
    logfile = "hippos.log"

    if len(sys.argv) > 1:
        config = sys.argv[1]
    else:
        print("HIPPOS config file not defined, using config.txt instead")
        print("To change the config file, use it as argument after HIPPOS")
        print("Example:")
        print("\t hippos <config file>\n")

    try:
        configread = open(config, 'r')
    except FileNotFoundError:
        print("The config file: '%s' can not be found" % (config))
        sys.exit(1)

    configlines = [line for line in configread]
    configread.close()

    for line in configlines:
        uncommented = line.split('#')[0]
        line_list = uncommented.split()

        if not line_list:
            continue
        option = line_list[0]

        if option == "docking_method":
            value = line_list[1]
            method = ["vina", "plants"]
            if value in method:
                docking_method = value
            else:
                print("docking method '%s' is not recognized" % (value))
                sys.exit(1)

        elif option == "docking_conf":
            docking_conf = line_list[1]

        elif option == "similarity_coef":
            value = line_list[1:]
            similarity_coef = value

        elif option == "simplified_ref":
            value = line_list[1:]
            simplified_ref = value

        elif option == "full_ref":
            value = line_list[1:]
            full_ref = value

        elif option == "full_nobb_ref":
            value = line_list[1:]
            full_nobb_ref = value

        elif option == "use_backbone":
            value = line_list[1]
            bool_val = ["yes", "no"]
            if value in bool_val:
                if value == "no":
                    use_backbone = False

        elif option == "residue_name":
            res_name = line_list[1:]
        elif option == "residue_number":
            res_num = line_list[1:]

        elif option == "res_weight1":
            res_weight1 = line_list[1:]
        elif option == "res_weight2":
            res_weight2 = line_list[1:]
        elif option == "res_weight3":
            res_weight3 = line_list[1:]
        elif option == "res_weight4":
            res_weight4 = line_list[1:]
        elif option == "res_weight5":
            res_weight5 = line_list[1:]

        elif option == "docking_score":
            value = line_list[1]
            bool_val = ["yes", "no"]
            if value in bool_val:
                if value == "no":
                    docking_score = False

        elif option == "output_mode":
            values = line_list[1:]
            mode = ["full", "full_nobb", "simplified"]
            for value in values:
                if value in mode:
                    output_mode[value] = True
                else:
                    print("output_mode '%s' is not recognized" % (value))
                    sys.exit(1)

        elif option == "simplified_outfile":
            simplified_outfile = line_list[1]
        elif option == "full_outfile":
            full_outfile = line_list[1]
        elif option == "full_nobb_outfile":
            full_nobb_outfile = line_list[1]
        elif option == "sim_outfile":
            sim_outfile = line_list[1]
        elif option == "logfile":
            logfile = line_list[1]
        elif option:
            print("Warning: '%s' option is not recognized" % (option))

    # Check output_mode value, if undefined then assign default value
    if output_mode == output_mode_undefined:
        output_mode['full'] = True

    parse_result = {
        'docking_method': docking_method,
        'docking_conf': docking_conf,
        'similarity_coef': similarity_coef,
        'full_ref': full_ref,
        'full_nobb_ref': full_nobb_ref,
        'simplified_ref': simplified_ref,
        'use_backbone': use_backbone,
        'residue_name': res_name,
        'residue_number': res_num,
        'res_weight1': res_weight1,
        'res_weight2': res_weight2,
        'res_weight3': res_weight3,
        'res_weight4': res_weight4,
        'res_weight5': res_weight5,
        'docking_score': docking_score,
        'output_mode': output_mode,
        'simplified_outfile': simplified_outfile,
        'full_outfile': full_outfile,
        'full_nobb_outfile': full_nobb_outfile,
        'sim_outfile': sim_outfile,
        'logfile': logfile,
    }
    '''import pprint
    pp = pprint.PrettyPrinter()
    pp.pprint(parse_result)'''
    return parse_result


def parse_config_genref():
    #  Default configuration
    config = "genref-config.txt"

    #  Identifiers initialization

    proteins = ""
    ligands = ""

    use_backbone = True
    res_name = []
    res_num = []
    res_weight1 = []
    res_weight2 = []
    res_weight3 = []
    res_weight4 = []
    res_weight5 = []

    output_mode = {"full": False, "full_nobb": False, "simplified": False}
    output_mode_undefined = {"full": False, "full_nobb": False, "simplified": False}

    outfile = "genref-results.txt"
    logfile = "hippos-genref.log"

    if len(sys.argv) > 1:
        config = sys.argv[1]
    else:
        print("HIPPOS-GENREF config file not defined, using genref-config.txt instead")
        print("To change the config file, use it as argument after HIPPOS")
        print("Example:")
        print("\t hippos-genref <config file>\n")

    try:
        configread = open(config, 'r')
    except FileNotFoundError:
        print("The config file: '%s' can not be found" % (config))
        sys.exit(1)

    configlines = [line for line in configread]
    configread.close()

    for line in configlines:
        uncommented = line.split('#')[0]
        line_list = uncommented.split()

        if not line_list:
            continue
        option = line_list[0]

        if option == "proteins":
            proteins = line_list[1:]
        elif option == "ligands":
            ligands = line_list[1:]

        elif option == "use_backbone":
            value = line_list[1]
            bool_val = ["yes", "no"]
            if value in bool_val:
                if value == "no":
                    use_backbone = False

        elif option == "residue_name":
            res_name = line_list[1:]
        elif option == "residue_number":
            res_num = line_list[1:]

        elif option == "res_weight1":
            res_weight1 = line_list[1:]
        elif option == "res_weight2":
            res_weight2 = line_list[1:]
        elif option == "res_weight3":
            res_weight3 = line_list[1:]
        elif option == "res_weight4":
            res_weight4 = line_list[1:]
        elif option == "res_weight5":
            res_weight5 = line_list[1:]

        elif option == "output_mode":
            values = line_list[1:]
            mode = ["full", "full_nobb", "simplified"]
            for value in values:
                if value in mode:
                    output_mode[value] = True
                else:
                    print("output_mode '%s' is not recognized" % (value))
                    sys.exit(1)

        elif option == "outfile":
            outfile = line_list[1]
        elif option == "logfile":
            logfile = line_list[1]
        elif option:
            print("Warning: '%s' option is not recognized" % (option))

    # Check output_mode value, if undefined then assign default value
    if output_mode == output_mode_undefined:
        output_mode['full'] = True

    parse_result = {
        'proteins': proteins,
        'ligands': ligands,
        'use_backbone': use_backbone,
        'residue_name': res_name,
        'residue_number': res_num,
        'res_weight1': res_weight1,
        'res_weight2': res_weight2,
        'res_weight3': res_weight3,
        'res_weight4': res_weight4,
        'res_weight5': res_weight5,
        'outfile': outfile,
        'logfile': logfile,
        'output_mode': output_mode,
    }

    return parse_result
