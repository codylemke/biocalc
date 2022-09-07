#! /usr/bin/python
"""

"""
import sys
import argparse
from Bio import AlignIO
from Bio import SeqIO
import numpy as np
from primer_design.optimize import genetic_algorithm
from pprint import pprint
from primer_design.codon_tables import CodonTable
from Bio.Seq import Seq

import pandas as pd


def find_seq_in_msa(seq, msa):
    """
      seq is a SeqRecord object
      msa is MultipleSequenceAlignment object
    """

    out = None
    for m in msa:
        if m.name == seq.name:
            return m
    return out


class PeptideCodonGapped:
    def __init__(self, name, nucleotide_seq, aa_seq, codon_table_id=1):
        self.name = name
        self.nucleotide_seq = nucleotide_seq
        self.aa_seq = aa_seq
        self.codon_table = CodonTable(codon_table_id)
        self.aligned_codons = list()
        nucl_pos = 0
        for aa in self.aa_seq:
            if aa != '-':
                codon = self.nucleotide_seq[nucl_pos:nucl_pos + 3]
                if self.codon_table[codon] == aa:
                    self.aligned_codons.append(codon)
                else:
                    raise ValueError("aa and nucl sequences do not correspond.")
                nucl_pos += 3
            else:
                self.aligned_codons.append(None)

    def __getitem__(self, index):
        return self.aligned_codons[index]

    def get_aa(self, index):
        return self.aa_seq[index]


def select_codons(gapped_frag, reference_optimizations, optimization_strategy, gapped=False, codon_table_id=1):
    out = ""
    codon_table = CodonTable(codon_table_id)
    for i, p in enumerate(gapped_frag):
        if p != '-':
            codon = None

            # check if the aa is in any reference optimization
            for ref_opt in reference_optimizations:
                if ref_opt[i] is not None:
                    if ref_opt.get_aa(i) == p:
                        codon = ref_opt[i]
                        break
            if codon is None:
                if optimization_strategy == "random":
                    codon_options = codon_table.aa_to_codons(p)
                    c_i = np.random.randint(len(codon_options))
                    codon = codon_options[c_i]
            out = out + codon
        else:  # p == '-'
            if gapped:
                out = out + '---'
    return out


def write_optimized_fasta(frag_df, output_handle, outer_left_overhang, outer_right_overhang, inner_left_overhang,
                          inner_right_overhang):
    if outer_left_overhang is None:
        outer_left_overhang = ""
    if outer_right_overhang is None:
        outer_right_overhang = ""

    # inner_left_overhang = ""
    # inner_right_overhang = ""
    # inner_right_overhang = str(Seq(inner_left_overhang, generic_dna).reverse_complement().seq)
    if inner_left_overhang is None:
        inner_left_overhang = ""
    if inner_right_overhang is None:
        inner_right_overhang = ""

    max_position = frag_df['position'].max()

    for i, r in frag_df.iterrows():
        out_seq = r['optimized_sequence']
        if r['position'] == 0:  # it's a fragment for single-fragment assemblies
            out_seq = outer_left_overhang + out_seq + outer_right_overhang
        if r['position'] == 1:  # it's the first part of a multi-fragment assembly
            out_seq = outer_left_overhang + out_seq
        if r['position'] == max_position:  # it's the last part of a multi-fragment assembly
            out_seq = out_seq + outer_right_overhang

        if r['position'] > 1:  # it's part of a multi-fragment assembly, but not the first part
            out_seq = inner_left_overhang + out_seq

        if ((r['position'] > 0) and (
                r['position'] < max_position)):  # it's part of a multi-fragment assembly, but not the last part
            out_seq = out_seq + outer_right_overhang

        print(">" + r.name + "\n" + out_seq, file=output_handle)


def main(input_handle, output_handle, params):
    codon_table = 1
    # read the fragments
    frag_df = pd.read_csv(input_handle, sep="\t", index_col=0)

    reference_optimizations = list()
    if params.reference_optimizations is not None:
        ref_opts = [x for x in SeqIO.parse(params.reference_optimizations, "fasta")]
        ref_alignment = AlignIO.read(params.reference_msa, params.reference_msa_type)

        for seq in ref_opts:
            aligned_aa = find_seq_in_msa(seq, ref_alignment)
            if aligned_aa is None:
                raise ValueError("Could not find alignment for reference optimization %s" % seq.name)
            reference_optimizations.append(PeptideCodonGapped(seq.name, str(seq.seq), str(aligned_aa.seq), codon_table))

    frag_df["optimized_sequence"] = frag_df['gapped sequence'].map(
        lambda gapped_frag: select_codons(gapped_frag, reference_optimizations, params.optimization_strategy,
                                          params.gapped_output, codon_table))

    write_optimized_fasta(frag_df, output_handle, params.outer_left_overhang, params.outer_right_overhang,
                          params.inner_left_overhang, params.inner_right_overhang)


####TODO: FIGURE OUT OVERHANGS!!!!
#### possibly with overhangs by position,

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # this should be the fragments file
    parser.add_argument("-i", default=None,
                        help="A file listing the fragments that need to be assembled, as amino-acid sequences."
                             "First column is name, second column is the number of times the fragment is used, third "
                             "column is position, fourth column is sequence, fifth column is gapped sequence. "
                             "If position = 0, that indicates that the fragment represents a whole sequence and is not "
                             "to be combined with other fragments, otherwise indicates fragment "
                             "order in multipart assemblies.")
    parser.add_argument("-o", default=None, help="Optimized fragment sequences will be written here, in fasta format.")
    # TODO: maybe encode position in fragment name?
    parser.add_argument("--gapped_output", action='store_true', default=False, required=False,
                        help="If true, the output will gapped according to the input.")
    parser.add_argument("--reference_msa", type=str, default=None, required=False,
                        help="A fasta protein alignment file.")  # TODO: maybe encode position in fragment name?
    parser.add_argument("--reference_msa_type", default="fasta", help="type of alignment file of the reference msa.")
    parser.add_argument("--reference_optimizations", type=str, required=False,
                        help="Fasta file containing codon optimized versions of some of the proteins in the alignment. "
                             "If amino acids in fragments match the amino acids in these sequences, then they will be "
                             "used for the output. Priority will go to the top sequence in this file and then in "
                             "order after that.")
    parser.add_argument("--outer_left_overhang", type=str, default=None, required=False,
                        help="A nucleotide sequence to add to the 5' end of outer fragments.")
    parser.add_argument("--outer_right_overhang", type=str, default=None, required=False,
                        help="A nucleotide sequence to add to the 3' end of outer fragments.")
    parser.add_argument("--inner_left_overhang", type=str, default=None, required=False,
                        help="A nucleotide sequence to add to the 5' end of inner fragments.")
    parser.add_argument("--inner_right_overhang", type=str, default=None, required=False,
                        help="A nucleotide sequence to add to the 3' end of internal fragments.")
    # parser.add_argument("--fragments", type=str, required=True,
    #                     help="A file listing the fragments that need to be assembled, as amino-acid sequences. "
    #                          "First column is name, second column is the number of times the fragment is used, third "
    #                          "column is position, fourth column is sequence, fifth column is gapped sequence. "
    #                          "If position = 0, that indicates that the fragment represents a whole sequence and is "
    #                          "not to be combined with other fragments, otherwise indicates fragment order "
    #                          "in multipart assemblies.") #TODO: maybe encode position in fragment name?
    parser.add_argument("--optimization_strategy", type=str, default="random", choices=['random'],
                        help="How to choose codons for positions where there is no reference codon.")

    # TODO: HINT if you want to enforce some overhang sqeuence, then make a reference optimization containing just
    #  the border positions and put it at the top of the reference_optimization file use the appropriate
    # TODO: I think by using the golden-hinges and DNA chisel libraries, we can make this program fairly versatile,
    #  including selection of overhang sequences and adding or removing restriction sites, etc.

    args = parser.parse_args()

    if args.i is not None:
        input_handle = open(args.i, "r")
    else:
        input_handle = sys.stdin

    if args.o is not None:
        output_handle = open(args.o, "w")
    else:
        output_handle = sys.stdout

    main(input_handle, output_handle, args)
