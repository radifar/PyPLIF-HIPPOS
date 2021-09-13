#!/usr/bin/env python

from __future__ import print_function

from time import time

from initialize.parse_conf import ParseConfigGenref
from ifp_processing import get_refbitstring


def main():
    x = time()
    genref_config = ParseConfigGenref()
    genref_config.parse_config()
    proteins = genref_config.proteins
    ligands = genref_config.ligands

    # calculate bitstring
    bitstrings = get_refbitstring(genref_config)

    # write bitstring result to output file
    outfile = open(genref_config.outfile, "w")

    # A nested for loop, the outer one evaluate protein ligand pair
    # the inner one evaluate every residue.
    for protein, ligand, bitstring in zip(proteins, ligands, bitstrings):
        # for every protein ligand pair initialize the bits
        full_bits = ""
        full_nobb_bits = ""
        simp_bits = ""

        # for every residue convert the bits to string then concatenate with bitstring
        for resname in genref_config.residue_name:
            if genref_config.output_mode["full"]:
                full_bits += bitstring[resname].full_bit.to01()
            if genref_config.output_mode["full_nobb"]:
                full_nobb_bits += bitstring[resname].nobb_bit.to01()
            if genref_config.output_mode["simplified"]:
                simp_bits += bitstring[resname].simp.to01()

        # augment the bits with description and new line
        if full_bits:
            full_bits = "full_bitstring: " + full_bits + "\n"
        if full_nobb_bits:
            full_nobb_bits = "full_nobb_bitstring: " + full_nobb_bits + "\n"
        if simp_bits:
            simp_bits = "simplified_bitstring: " + simp_bits + "\n"

        # dump everything to output file
        outfile.write(
            "Protein: %s\nLigand:  %s\n%s%s%s\n"
            % (protein, ligand, full_bits, full_nobb_bits, simp_bits)
        )

    y = time()
    total_time = "\nTotal time taken %.3f s." % (y - x)
    outfile.write(total_time)
    outfile.close()

    print(total_time)


if __name__ == "__main__":
    main()
