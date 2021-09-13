from bitarray import bitarray
from SIMILARITY_FORMULA import sim_dict

"""
Based on "Similarity Coefficients for Binary Chemoinformatics
Data: Overview and Extended Comparison Using Simulated and 
Real Data Sets" by Todeschini et al. 2012
dx.doi.org/10.1021/ci300261r
"""


def clean_omitted_interactions(refbit, tgtbit):
    clean_refbit = ""
    clean_tgtbit = ""
    for ref, tgt in zip(refbit, tgtbit):
        if tgt != "n":
            clean_refbit += ref
            clean_tgtbit += tgt
    return clean_refbit, clean_tgtbit


def count_abcdp(refbit, tgtbit):
    refbit, tgtbit = clean_omitted_interactions(refbit, tgtbit)
    refbit = bitarray(refbit)
    tgtbit = bitarray(tgtbit)
    a = (refbit & tgtbit).count()
    b = (~refbit & tgtbit).count()
    c = (refbit & ~tgtbit).count()
    d = (~refbit & ~tgtbit).count()
    p = a + b + c + d

    return a, b, c, d, p


def how_similar(abcdp, sim_coef):
    a, b, c, d, p = abcdp  # lgtm [py/unused-local-variable]
    try:
        sim = eval(sim_dict[sim_coef])
    except ZeroDivisionError:
        sim = "NA"
    return sim


def replace_bit_char(bitstring, bit_index_list):
    for i, v in enumerate(bit_index_list):
        if v == 1:
            bitstring = bitstring[:i] + "n" + bitstring[i+1:]
    return bitstring
