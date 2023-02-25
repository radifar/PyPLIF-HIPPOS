#!/usr/bin/env python

from __future__ import print_function

import sys
from time import time
from typing import Type, List, Tuple, Dict

from initialize.parse_conf import ParseConfig
from initialize.parse_docking_conf import parse_plants_conf, parse_vina_conf
from parse_mol import (
    parse_complex,
    parse_ligands,
    parse_protein,
    enumerate_ligand_files,
    enumerate_ligand_file_list,
)
from ifp_processing import (
    Residue,
    get_bitstring,
    get_direct_bitstring,
    get_complex_bitstring,
)
from similarity import count_abcdp, how_similar, replace_bit_char
from observer import setup_dict, do_task


def collect_ligand(ligand_pose: List[str], config: Type[ParseConfig]) -> None:
    """Process file ligand names and add it to ligand_pose

    Parameters
    ----------
    ligand_pose : List[str]
        list consisting ligand file names
    config : ParseConfig
        ParseConfig object which contain configuration data
    """

    ligand_files = config.ligand_files
    multiple_ligand_list = config.multiple_ligand_list
    if ligand_files:
        enumerate_ligand_files(ligand_pose, ligand_files)
    if multiple_ligand_list:
        enumerate_ligand_file_list(ligand_pose, multiple_ligand_list)


def process_input_to_bitstring(
    config: Type[ParseConfig],
) -> Tuple[List[str], List[str], Dict[str, Type[Residue]]]:
    """Using provided configuration generate ligand_pose, scorelist, and bitstrings

    Parameters
    ----------
    config : ParseConfig
        ParseConfig object which contain configuration data

    Returns
    -------
    ligand_pose, scorelist, bitstrings: Tuple[List[str], List[str], Dict[str, Type[Residue]]]
        Generate ligand_pose, scorelist, and bitstrings
    """

    ligand_pose = []
    scorelist = []

    if config.direct_ifp:
        if not config.complex_list:
            collect_ligand(ligand_pose, config)
            print(ligand_pose)
            ligand_mol_list = parse_ligands(ligand_pose)
            protein_mol = parse_protein(config.protein)

            bitstrings = get_direct_bitstring(protein_mol, ligand_mol_list, config)
            scorelist = [""] * len(ligand_pose)
        else:
            complex_list = []
            with open(config.complex_list, "r") as complex:
                for pair in complex:
                    complex_list.append(pair)
            scorelist = [""] * len(complex_list)
            ligand_pose = complex_list

            complex_mol_list = parse_complex(complex_list)
            bitstrings = get_complex_bitstring(complex_mol_list, config)
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

        do_task("is_docking_output_missing", docking_results)

        scorelist = docking_results["scorelist"]
        bitstrings = get_bitstring(docking_results, config)

    return ligand_pose, scorelist, bitstrings


def main():
    x = time()
    config = ParseConfig()
    config.parse_config()

    for mode in config.output_mode:
        if config.output_mode[mode]:
            setup_dict[mode](config)

    setup_dict["logfile"](config)
    if config.similarity_coef:
        setup_dict["similarity"](config)

    ligand_pose, scorelist, bitstrings = process_input_to_bitstring(config)

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
            residue = bitstrings[resname]
            bit_replace_index = residue.bit_replace_index
            simp_bit_replace_index = residue.simp_bit_replace_index

            if config.output_mode["simplified"]:
                simp_res_bit = residue.simp_bits_list[pose].to01()
                if bool(sum(simp_bit_replace_index)):
                    simp_res_bit = replace_bit_char(
                        simp_res_bit, simp_bit_replace_index
                    )
                simp_bits += simp_res_bit
                if pose == 0:
                    bitlength = len(simp_res_bit)
                    bit_end = bit_start + bitlength - 1
                    do_task(
                        "simplified_bit_log", resname, bitlength, bit_start, bit_end
                    )
                    bit_start += bitlength

            if config.output_mode["full"]:
                full_res_bit = residue.full_bits_list[pose].to01()
                if bool(sum(bit_replace_index)):
                    full_res_bit = replace_bit_char(full_res_bit, bit_replace_index)
                full_bits += full_res_bit

            if config.output_mode["full_nobb"]:
                nobb_res_bit = residue.full_nobb_list[pose].to01()
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
