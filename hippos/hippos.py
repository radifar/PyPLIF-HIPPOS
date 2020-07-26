#!/usr/bin/env python

from __future__ import print_function

import sys
from time import time

from initialize.parse_conf import parse_config
from ifp_processing import get_bitstring
from similarity import count_abcdp, how_similar


def main():
    x = time()

    hippos_config = parse_config()

    logfile = open(hippos_config['logfile'], 'w')  # Output #4
    '''
    Parse docking configuration file
    Get docking results:
    protein			==> OBMol Object
    docked_ligands 	==> List of OBMol
    docked_proteins ==> List of OBMol (only for PLANTS)
    mollist			==> List of ligand name + pose number
    scorelist		==> List of docking score
    '''

    if hippos_config['docking_method'] == "plants":
        from initialize.parse_docking_conf import parse_plants_conf
        docking_conf = hippos_config['docking_conf']
        docking_results = parse_plants_conf(docking_conf)
    else:
        from initialize.parse_docking_conf import parse_vina_conf
        docking_conf = hippos_config['docking_conf']
        docking_results = parse_vina_conf(docking_conf)

    if len(docking_results['docked_ligands']) == 0:
        print("The docking output could not be found. Please check your docking result.")
        logfile.write("The docking output could not be found. Please check your docking result.")
        logfile.close()
        sys.exit(1)
    '''
    Get Bitstring using docking results & hippos configuration
    bitstrings		==> Dictionary, with resname as key
                        Residue object as value
    '''

    bitstrings = get_bitstring(docking_results, hippos_config)
    '''
    Write Output & Log files
    '''

    scorelist = docking_results['scorelist']
    ligand_pose = []
    if hippos_config['docking_method'] == "plants":
        for mol in docking_results['mollist']:
            mol = mol.split('_')
            new_name = mol[0] + '_' + mol[-1]
            ligand_pose.append(new_name)
    if hippos_config['docking_method'] == "vina":
        ligand_pose = docking_results['mollist']

    output_mode = hippos_config['output_mode']
    simplified_flag = output_mode['simplified']
    full_flag = output_mode['full']
    full_nobb_flag = output_mode['full_nobb']

    if simplified_flag:
        simplified_outfile = open(hippos_config['simplified_outfile'], 'w')  # Output #1
    if full_flag:
        full_outfile = open(hippos_config['full_outfile'], 'w')  # Output #2
    if full_nobb_flag:
        full_nobb_outfile = open(hippos_config['full_nobb_outfile'], 'w')  # Output #3


    logfile.write('Ligand name is %s with %s poses\n\n' % \
        (ligand_pose[0].split('_')[0], len(ligand_pose))) # Output Logfile

    similarity_coef = hippos_config['similarity_coef']
    if similarity_coef:
        sim_outfile = open(hippos_config['sim_outfile'], 'w')  # Output #5
        logfile.write('similarity coefficient used are %s\n' % \
            (', '.join(similarity_coef)))     # Output Logfile
    if simplified_flag:
        logfile.write('%s %s %s %s\n' % \
            ('RESNAME','length','startbit','endbit'))     # Output Logfile
    '''
    Iterate through pose and write the ligand+pose, score, bitstring
    '''

    pose = 0
    log_flag = True
    bitstring_zero = False

    for ligand_name, score in zip(ligand_pose, scorelist):
        ligand_name = ligand_name.replace(' ', '_')
        simp_bits = ''
        full_bits = ''
        nobb_bits = ''

        bit_start = 1
        for resname in hippos_config['residue_name']:
            if simplified_flag:
                simp_res_bit = bitstrings[resname].simp_bits_list[pose].to01()
                simp_bits += simp_res_bit
            if full_flag:
                full_res_bit = bitstrings[resname].full_bits_list[pose].to01()
                full_bits += full_res_bit
            if full_nobb_flag:
                nobb_res_bit = bitstrings[resname].full_nobb_list[pose].to01()
                nobb_bits += nobb_res_bit
            if log_flag & simplified_flag:
                bitlength = len(simp_res_bit)
                bit_end = bit_start + bitlength - 1
                logfile.write('%-10s %-6s %-7s %s\n' % \
                    (resname, bitlength, bit_start, bit_end)) # Output Logfile
                bit_start += bitlength
        log_flag = False

        if simplified_flag:
            simplified_outfile.write('%s %s %s\n' % \
                (ligand_name.ljust(16), score.ljust(9), simp_bits)) # Output Simplified
        if full_flag:
            full_outfile.write('%s %s %s\n' % \
                (ligand_name.ljust(16), score.ljust(9), full_bits)) # Output Full
        if full_nobb_flag:
            full_nobb_outfile.write('%s %s %s\n' % \
                (ligand_name.ljust(16), score.ljust(9), nobb_bits)) # Output No BB
        '''
        If similarity coef requested => calculate abcd and p
        '''

        if similarity_coef:
            abcdp_list = []
            coefficient = []
            if full_flag:
                for full in hippos_config['full_ref']:
                    abcdp_list.append(count_abcdp(full, full_bits))
            elif full_nobb_flag:
                for nobb in hippos_config['full_nobb_ref']:
                    abcdp_list.append(count_abcdp(nobb, nobb_bits))
            else:
                for simp in hippos_config['simplified_ref']:
                    abcdp_list.append(count_abcdp(simp, simp_bits))
            for sim_coef in similarity_coef:
                for abcdp in abcdp_list:
                    similarity_value = how_similar(abcdp, sim_coef)
                    try:
                        coefficient.append('%.3f' % similarity_value)
                    except TypeError:
                        coefficient.append('%s' % similarity_value)
                        bitstring_zero = True
            sim_outfile.write('%s %s\n' % \
                (ligand_name.ljust(16), ' '.join(coefficient))) # Output Similarity

        pose += 1
    '''
    Close all file
    '''

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
        print("It appears that one of the target or reference bitstring is zero, \
the similarity value will be set to 'NA'. Check the ligand pose that generate 'NA' value.\n")
        logfile.write('''
It appears that one of the target or reference bitstring is zero,
Check the ligand pose that generate 'NA' value.
        ''')
    print('Total time taken %.3f s.' % z)
    logfile.write('\nTotal time taken %.3f s.' % z)
    logfile.close()


if __name__ == "__main__":
    main()
