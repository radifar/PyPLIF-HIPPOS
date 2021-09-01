#!/usr/bin/env python

from __future__ import print_function

import sys
from time import time

from initialize.parse_conf import ParseConfig
from initialize.parse_docking_conf import parse_plants_conf, parse_vina_conf
from parse_mol import (
    parse_ligands,
    parse_protein,
    enumerate_ligand_files,
    enumerate_ligand_file_list
)
from ifp_processing import (
    get_bitstring,
    get_direct_bitstring,
    get_complex_bitstring
)
from similarity import count_abcdp, how_similar, replace_bit_char


def main():
    x = time()
    hippos_config = ParseConfig()
    hippos_config.parse_config()

    ligand_pose = []
    scorelist = []

    if hippos_config.direct_ifp:
        if not hippos_config.complex_list:
            ligand_files = hippos_config.ligand_files
            ligand_file_list = hippos_config.ligand_file_list
            if ligand_files:
                enumerate_ligand_files(ligand_pose, ligand_files)
            if ligand_file_list:
                enumerate_ligand_file_list(ligand_pose, ligand_file_list)

            ligand_mol_list = parse_ligands(ligand_pose)
            protein_mol = parse_protein(hippos_config.protein)

            bitstrings = get_direct_bitstring(protein_mol, ligand_mol_list, hippos_config)
            for pose in ligand_pose:
                scorelist.append("")
        else:
            bitstrings = get_complex_bitstring("complex_list", hippos_config)
            print(bitstrings)
    else:
        if hippos_config.docking_method == "plants":
            docking_results = parse_plants_conf(hippos_config.docking_conf)
            for mol in docking_results["ligand_pose"]:
                mol = mol.split("_")
                new_name = mol[0] + "_" + mol[-1]
                ligand_pose.append(new_name)
        else:
            docking_results = parse_vina_conf(hippos_config.docking_conf)
            ligand_pose = docking_results["ligand_pose"]

        # checking docking output, if not found then exit.
        if len(docking_results["docked_ligands"]) == 0:
            missing_docking_output = (
                "The docking output could not be found. Please check your docking result."
            )

            print(missing_docking_output)
            with open(hippos_config.logfile, "w") as logfile:
                logfile.write(missing_docking_output)

            sys.exit(1)

        scorelist = docking_results["scorelist"]
        bitstrings = get_bitstring(docking_results, hippos_config)

    # set flag for every chosen output mode
    output_mode = hippos_config.output_mode
    simplified_flag = output_mode["simplified"]
    full_flag = output_mode["full"]
    full_nobb_flag = output_mode["full_nobb"]

    # set file handler for every chosen output mode
    if simplified_flag:
        simplified_outfile = open(hippos_config.simplified_outfile, "w")  # Output #1
    if full_flag:
        full_outfile = open(hippos_config.full_outfile, "w")  # Output #2
    if full_nobb_flag:
        full_nobb_outfile = open(hippos_config.full_nobb_outfile, "w")  # Output #3

    # write ligand info and similarity coef info
    logfile = open(hippos_config.logfile, "w")
    logfile.write(
        "Ligand name is %s with %s poses\n\n"
        % (ligand_pose[0].split("_")[0], len(ligand_pose))
    )  # Output Logfile

    similarity_coef = hippos_config.similarity_coef
    if similarity_coef:
        sim_outfile = open(hippos_config.sim_outfile, "w")  # Output #5
        logfile.write(
            "similarity coefficient used are %s\n" % (", ".join(similarity_coef))
        )  # Output Logfile

    # if simplified then write the length and position for each bitstring
    if simplified_flag:
        logfile.write(
            "%s %s %s %s\n" % ("RESNAME", "length", "startbit", "endbit")
        )  # Output Logfile

    # Iterate through pose and write the ligand+pose, docking score,
    # similarity coef, bitstring
    log_flag = True
    bitstring_zero = False

    for pose, (ligand_name, score) in enumerate(zip(ligand_pose, scorelist)):
        ligand_name = ligand_name.replace(" ", "_").ljust(16)
        score = score.ljust(9)
        simp_bits = ""
        full_bits = ""
        nobb_bits = ""

        # Concatenate bitstring from every residue, then write to their respective
        # output file
        bit_start = 1
        for resname in hippos_config.residue_name:
            bit_replace_index = bitstrings[resname].bit_replace_index
            simp_bit_replace_index = bitstrings[resname].simp_bit_replace_index
            if simplified_flag:
                simp_res_bit = bitstrings[resname].simp_bits_list[pose].to01()
                if bool(sum(simp_bit_replace_index)):
                    simp_res_bit = replace_bit_char(simp_res_bit, simp_bit_replace_index)
                simp_bits += simp_res_bit
            if full_flag:
                full_res_bit = bitstrings[resname].full_bits_list[pose].to01()
                if bool(sum(bit_replace_index)):
                    full_res_bit = replace_bit_char(full_res_bit, bit_replace_index)
                full_bits += full_res_bit
            if full_nobb_flag:
                nobb_res_bit = bitstrings[resname].full_nobb_list[pose].to01()
                if bool(sum(bit_replace_index)):
                    nobb_res_bit = replace_bit_char(nobb_res_bit, bit_replace_index)
                nobb_bits += nobb_res_bit
            if log_flag & simplified_flag:
                bitlength = len(simp_res_bit)
                bit_end = bit_start + bitlength - 1
                logfile.write(
                    "%-10s %-6s %-7s %s\n" % (resname, bitlength, bit_start, bit_end)
                )  # Output Logfile
                bit_start += bitlength
        log_flag = False

        if simplified_flag:
            simplified_outfile.write("%s %s %s\n" % (ligand_name, score, simp_bits))
        if full_flag:
            full_outfile.write("%s %s %s\n" % (ligand_name, score, full_bits))
        if full_nobb_flag:
            full_nobb_outfile.write("%s %s %s\n" % (ligand_name, score, nobb_bits))

        # If similarity coef requested => calculate abcd and p
        if similarity_coef:
            abcdp_list = []
            coefficient = []
            if full_flag:
                for full in hippos_config.full_ref:
                    abcdp_list.append(count_abcdp(full, full_bits))
            elif full_nobb_flag:
                for nobb in hippos_config.full_nobb_ref:
                    abcdp_list.append(count_abcdp(nobb, nobb_bits))
            else:
                for simp in hippos_config.simplified_ref:
                    abcdp_list.append(count_abcdp(simp, simp_bits))
            for sim_coef in similarity_coef:
                for abcdp in abcdp_list:
                    similarity_value = how_similar(abcdp, sim_coef)
                    try:
                        coefficient.append("%.3f" % similarity_value)
                    except TypeError:
                        coefficient.append("%s" % similarity_value)
                        bitstring_zero = True
            sim_outfile.write(
                "%s %s\n" % (ligand_name, " ".join(coefficient))
            )  # Output Similarity

    # Close all file
    if simplified_flag:
        simplified_outfile.close()
    if full_flag:
        full_outfile.close()
    if full_nobb_flag:
        full_nobb_outfile.close()
    if similarity_coef:
        sim_outfile.close()

    y = time()
    z = y - x
    if bitstring_zero:
        bitstring_error = """
It appears that one of the target or reference bitstring is zero,
Check the ligand pose that generate 'NA' value.
        """
        print(bitstring_error)
        logfile.write(bitstring_error)
    print("Total time taken %.3f s." % z)
    logfile.write("\nTotal time taken %.3f s." % z)
    logfile.close()


if __name__ == "__main__":
    main()
