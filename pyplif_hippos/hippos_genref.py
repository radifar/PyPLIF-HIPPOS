#!/usr/bin/env python

from __future__ import print_function

from time import time

from initialize.parse_conf import parse_config_genref
from ifp_processing import get_refbitstring

def main():
    x = time()
    genref_config = parse_config_genref()
    proteins = genref_config['proteins']
    ligands = genref_config['ligands']

    # calculate bitstring
    bitstrings = get_refbitstring(genref_config)

    # write bitstring result to output file
    outfile = open(genref_config['outfile'], 'w')
    
    for protein, ligand, bitstring in zip(proteins, ligands, bitstrings):
        full_bits = ''
        full_nobb_bits = ''
        simp_bits = ''

        for resname in genref_config['residue_name']:
            if genref_config['output_mode']['full']:
                full_bits += bitstring[resname].full_bit.to01()
            if genref_config['output_mode']['full_nobb']:
                full_nobb_bits += bitstring[resname].nobb_bit.to01()
            if genref_config['output_mode']['simplified']:
                simp_bits += bitstring[resname].simp.to01()
        if full_bits:
            full_bits = 'full_bitstring: ' + full_bits + '\n'
        if full_nobb_bits:
            full_nobb_bits = 'full_nobb_bitstring: ' + full_nobb_bits + '\n'
        if simp_bits:
            simp_bits = 'simplified_bitstring: ' + simp_bits + '\n'
        outfile.write('Protein: %s\nLigand:  %s\n%s%s%s\n' % (protein, ligand, full_bits, full_nobb_bits, simp_bits))

    y = time()
    total_time = '\nTotal time taken %.3f s.' % (y - x)
    outfile.write(total_time)
    outfile.close()

    print(total_time)


if __name__ == "__main__":
    main()
