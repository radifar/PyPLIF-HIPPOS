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
from observer import setup_dict, do_task


def collect_ligand(ligand_pose, config):
    ligand_files = config.ligand_files
    ligand_file_list = config.ligand_file_list
    if ligand_files:
        enumerate_ligand_files(ligand_pose, ligand_files)
    if ligand_file_list:
        enumerate_ligand_file_list(ligand_pose, ligand_file_list)


def is_docking_output_missing(docking_results):
    if len(docking_results["docked_ligands"]) == 0:
        missing_docking_output = (
            "The docking output could not be found. Please check your docking result."
        )

        print(missing_docking_output)
        with open(config.logfile, "w") as logfile:
            logfile.write(missing_docking_output)

        sys.exit(1)


def process_input(config):
    ligand_pose = []
    scorelist = []

    if config.direct_ifp:
        if not config.complex_list:
            collect_ligand(ligand_pose, config)
            ligand_mol_list = parse_ligands(ligand_pose)
            protein_mol = parse_protein(config.protein)

            bitstrings = get_direct_bitstring(protein_mol, ligand_mol_list, config)
            scorelist = [""] * len(ligand_pose)
        else:
            bitstrings = get_complex_bitstring("complex_list", config)
            print(bitstrings)
    else:
        if config.docking_method == "plants":
            docking_results = parse_plants_conf(config.docking_conf)
            for mol in docking_results["ligand_pose"]:
                mol = mol.split("_")
                new_name = mol[0] + "_" + mol[-1]
                ligand_pose.append(new_name)
        else:
            docking_results = parse_vina_conf(config.docking_conf)
            ligand_pose = docking_results["ligand_pose"]

        is_docking_output_missing(docking_results)

        scorelist = docking_results["scorelist"]
        bitstrings = get_bitstring(docking_results, config)

    return ligand_pose, scorelist, bitstrings


def main():
    x = time()
    config = ParseConfig()
    config.parse_config()
    ligand_pose, scorelist, bitstrings = process_input(config)

    for mode in config.output_mode:
        if config.output_mode[mode]:
            setup_dict[mode](config)

    setup_dict["logfile"](config)
    if config.similarity_coef:
        setup_dict["similarity"](config)

    do_task("initial_information", config, ligand_pose)

    # Iterate through pose and write the ligand+pose, docking score,
    # similarity coef, bitstring
    for pose, (ligand_name, score) in enumerate(zip(ligand_pose, scorelist)):
        ligand_name = ligand_name.replace(" ", "_").ljust(16)
        score = score.ljust(9)
        simp_bits = full_bits = nobb_bits = ""

        # Concatenate bitstring from every residue, then write to their respective
        # output file
        bit_start = 1
        for resname in config.residue_name:
            bit_replace_index = bitstrings[resname].bit_replace_index
            simp_bit_replace_index = bitstrings[resname].simp_bit_replace_index
            if config.output_mode["simplified"]:
                simp_res_bit = bitstrings[resname].simp_bits_list[pose].to01()
                if bool(sum(simp_bit_replace_index)):
                    simp_res_bit = replace_bit_char(simp_res_bit, simp_bit_replace_index)
                simp_bits += simp_res_bit
                if pose == 0:
                    bitlength = len(simp_res_bit)
                    bit_end = bit_start + bitlength - 1
                    do_task("simplified_bit_log", resname, bitlength, bit_start, bit_end)
                    bit_start += bitlength

            if config.output_mode["full"]:
                full_res_bit = bitstrings[resname].full_bits_list[pose].to01()
                if bool(sum(bit_replace_index)):
                    full_res_bit = replace_bit_char(full_res_bit, bit_replace_index)
                full_bits += full_res_bit

            if config.output_mode["full_nobb"]:
                nobb_res_bit = bitstrings[resname].full_nobb_list[pose].to01()
                if bool(sum(bit_replace_index)):
                    nobb_res_bit = replace_bit_char(nobb_res_bit, bit_replace_index)
                nobb_bits += nobb_res_bit

        all_bits = (simp_bits, full_bits, nobb_bits)
        do_task("write_bitstrings", ligand_name, score, all_bits)

        # If similarity coef requested => calculate abcd and p
        if config.similarity_coef:
            abcdp_list = []
            coefficient = []
            if config.output_mode["full"]:
                for full in config.full_ref:
                    abcdp_list.append(count_abcdp(full, full_bits))
            elif config.output_mode["full_nobb"]:
                for nobb in config.full_nobb_ref:
                    abcdp_list.append(count_abcdp(nobb, nobb_bits))
            else:
                for simp in config.simplified_ref:
                    abcdp_list.append(count_abcdp(simp, simp_bits))
            for sim_coef in config.similarity_coef:
                for abcdp in abcdp_list:
                    similarity_value = how_similar(abcdp, sim_coef)
                    try:
                        coefficient.append("%.3f" % similarity_value)
                    except TypeError:
                        coefficient.append("%s" % similarity_value)
                        bitstring_error = (
                            "It appears that one of the target or reference bitstring is zero,\n"
                            "Check the ligand pose that generate 'NA' value."
                        )
                        do_task("bitstring_error", bitstring_error)
            do_task("write_similarity", ligand_name, coefficient)

    y = time()
    total_time = y - x
    print("Total time taken %.3f s." % total_time)
    do_task("write_time", total_time)
    do_task("close_files")


if __name__ == "__main__":
    main()
