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

    bitstrings = get_refbitstring(genref_config)

    outfile = open(genref_config['outfile'], 'w')
    '''
    if genref_config['use_backbone']:
      for protein, ligand, bitstring in zip(proteins, ligands, bitstrings):
          simp_bits = ''
          full_bits  = ''
          for resname in genref_config['residue_name']:
              simp_bits += bitstring[resname].simp.to01()
              full_bits  += bitstring[resname].full.to01()
              bb_res_bit = bitstring[resname].bb.to01()
              simp_bits += bb_res_bit
              full_bits  += bb_res_bit
          outfile.write('%s %s\n%s\n%s\n\n' % (protein, ligand, simp_bits, full_bits))
    else:
      for protein, ligand, bitstring in zip(proteins, ligands, bitstrings):
          simp_bits = ''
          full_bits  = ''
          for resname in genref_config['residue_name']:
              simp_bits += bitstring[resname].simp.to01()
              full_bits  += bitstring[resname].full.to01()
          outfile.write('%s %s\n%s\n%s\n\n' % (protein, ligand, simp_bits, full_bits))
    '''
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
    outfile.write('\nTotal time taken %.3f s.' % (y - x))
    outfile.close()

    print('Total time taken %.3f s.' % (y - x))


if __name__ == "__main__":
    main()
