#! /usr/bin/python
"""
  generates mutagenic primers for single codon or 2-consecutive mutations using, a Gibson assembly strategy.
  Input: A nucleotide sequence and a list of desired mutations.

  Enables primers to be made according to the following specification (and similar)
  F: 5'-(NNNNNNN)-(***)-(20-25N)-3'
  R: 3'-(20-25N)-(***)-(NNNNNNN)-5'
"""

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
import python_codon_tables as pct
from primer_design import parsers
import re

from dnachisel import DnaOptimizationProblem, AvoidPattern, EnzymeSitePattern, DnaNotationPattern, EnforceSequence, \
    EnforceTranslation, CodonOptimize

MIN_PRIMER_LENGTH = 15
MAX_PRIMER_LENGTH = 60


# INTERNAL_PRIMER_EXTENSION=10

def parse_target_mutations(mutations_string):
    parts = re.split(r'\s|_', mutations_string)
    out = list()
    if len(parts) > 0:
        for part in parts:
            match = re.match('^([A-Z])([0-9]+)([A-Z])$', part)
            if match:
                out.append({'position': int(match.group(2)), 'sourceaa': match.group(1).upper(),
                            'targetaa': match.group(3).upper()})
            else:
                raise ValueError("Error, cannot parse mutation specification %s" % part)
    else:
        raise ValueError("Error, cannot parse mutation specification %s" % mutations_string)
    out.sort(key=lambda x: x['position'])
    return out


def invert_codon_table(table):
    out = dict()
    for (aa, codons) in table.items():
        for codon in codons:
            out[codon] = aa
    return out


def sort_codon_table(table):
    out = dict()
    for (aa, codons) in table.items():
        out[aa] = [(codon, freq) for (codon, freq) in codons.items()]
        out[aa].sort(key=lambda x: x[1], reverse=True)
    return out


# TODO: be able to handle sequences where the mutations are close to the edge of the sequence,
#  and we therefore need to have the mutagenic primers cross into the vector sequence.
# TODO: might want to check for ambiguous binding.
def protein_position_to_codon(nucleotide_sequence, protein_position):
    """
    returns(position_of_first_base, codon)
    """
    pos = protein_position_to_nucleotide_position(protein_position)
    if pos + 3 > len(nucleotide_sequence):
        raise ValueError("Error, codon %d is outside the coding sequence" % protein_position)
    codon = nucleotide_sequence[pos:pos + 3]
    return pos, codon


def protein_position_to_nucleotide_position(protein_position):
    return (protein_position - 1) * 3


def find_primer(seq, target_tm, min_tm, max_tm, min_len, max_len, Tris, Mg, K, Na, dNTPs, primer_conc, gc_clamp):
    """
    finds a subsequence starting at the beginning of seq, with properties as specified by the rest of the parameters.
    if no sequence matching the constraints can be found, None is returned.
    """
    best_primer = dict()
    for l in range(min_len, max_len + 1):
        sub = seq[0:l]
        tm = MeltingTemp.Tm_NN(sub, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=6, dnac1=primer_conc / 2,
                               dnac2=primer_conc / 2)
        if (sub[-1] in "GC") or (not gc_clamp):
            if (tm >= min_tm) and (tm <= max_tm):
                if "seq" not in best_primer:
                    best_primer['seq'] = sub
                    best_primer['tm'] = tm
                elif abs(target_tm - tm) < abs(target_tm - best_primer['tm']):
                    best_primer['seq'] = sub
                    best_primer['tm'] = tm
    if len(best_primer) == 0:

        return None
    else:
        return best_primer


def main(inpt, output, params):
    target_mutations_file = params.target_mutations
    seqs = list(SeqIO.parse(inpt, "fasta"))
    if len(seqs) != 1:
        raise ValueError("Error: Input fasta must contain only one sequence")
    seq_name = seqs[0].name
    seq = seqs[0].seq

    if (len(seq) % 3) != 0:
        raise ValueError("Error: Input sequence must have length divisible by 3")

    mutations_file_lines = [x for x in parsers.SimpleParser(target_mutations_file)]
    target_mutant_names = dict()
    target_mutants = dict()
    for line in mutations_file_lines:
        parts = [x.strip() for x in line.split("\t")]

        if len(parts) == 1:
            target_mutants[parts[0]] = parse_target_mutations(parts[0])
            target_mutant_names[parts[0]] = seq_name + "_" + parts[0]
        elif len(parts) == 2:
            target_mutants[parts[1]] = parse_target_mutations(parts[1])
            target_mutant_names[parts[1]] = parts[0]
        else:
            print("Warning, too many fields in design file line: %s" % line, file=sys.stderr)
    #

    # dict where key is mutation specification, and value is the parsed mutation specification
    # (a lists [mutation number] of dicts {sourceaa, targetaa, position} )
    # target_mutants = {x: parse_target_mutations(x) for x in parsers.simple_parser(target_mutations_file)}

    new_primers = list()  # values are tuples of name, sequence
    # new_primer_names = dict()
    # assembly = dict()
    new_primer_tms = dict()

    # returns a dict where keys are AA and values are dicts where keys are codons and values are normalized frequency.
    aa_to_codons = pct.get_codons_table(
        params.target_species_taxid)
    if len(aa_to_codons) == 0:
        raise ValueError("Error: cannot find codon table corresponding to taxid:" % params.target_species_taxid)
    codon_to_aa = invert_codon_table(aa_to_codons)

    sorted_codon_table = sort_codon_table(aa_to_codons)

    # counters = {'left': 0, 'right': 0}

    # mutation_space=list() #for codon optimization

    for (mutant, mutations) in target_mutants.items():
        mutations.sort(key=lambda x: x['position'])
        left_edge = protein_position_to_nucleotide_position(mutations[0]['position'])
        right_edge = protein_position_to_nucleotide_position(mutations[-1]['position']) + 3
        # insert_sequence = seq[left_edge:right_edge]
        modified_sequence = str(seq)
        for i in range(len(mutations)):
            (nuc_pos, codon) = protein_position_to_codon(str(seq), mutations[i]["position"])
            codon_aa = codon_to_aa[codon]
            if codon_aa != mutations[i]["sourceaa"]:
                print("Warning: mutation %s does not match template DNA %s -> %s" % (
                  mutations[i]['sourceaa'] + str(mutations[i]["position"]), codon, codon_aa), file=sys.stderr)
            new_codon = sorted_codon_table[mutations[i]['targetaa']][0][0]
            modified_sequence = modified_sequence[0:nuc_pos] + new_codon + modified_sequence[nuc_pos + 3:]
        modified_sequence = Seq(modified_sequence)

        # TODO: right here is where we'd check modified sequence to make
        #  sure it doesn't have any forbidden restriction sites.
        if (params.enzymes_to_avoid is not None) or (params.motifs_to_avoid is not None):
            constraints = list()
            constraints.append(EnforceSequence(sequence=modified_sequence[0:left_edge], location=(0, left_edge, 1)))
            constraints.append(EnforceSequence(sequence=modified_sequence[right_edge:],
                                               location=(right_edge, len(modified_sequence), 1)))
            constraints.append(EnforceTranslation(genetic_table='Standard', location=(0, len(modified_sequence), 1),
                                                  start_codon='keep'))
            # TODO: Handle if genetic table is not Standard, also handle if this is a
            #  subsequence of a larger sequence, possibly handled via annotations.
            if params.enzymes_to_avoid is not None:
                for s in params.enzymes_to_avoid.split(","):
                    s = s.strip()
                    constraints.append(AvoidPattern(EnzymeSitePattern(s)))
            if params.motifs_to_avoid is not None:
                for s in params.motifs_to_avoid.split(","):
                    s = s.strip()
                    constraints.append(AvoidPattern(DnaNotationPattern(s)))
            objectives = [CodonOptimize(species=None, method='use_best_codon', location=(0, len(modified_sequence), 1),
                                        codon_usage_table=aa_to_codons, boost=1.0)]

            opt_prob = DnaOptimizationProblem(str(modified_sequence), constraints=constraints, objectives=objectives,
                                              logger=None, mutation_space=None)
            if not opt_prob.all_constraints_pass():
                print(mutant)
                print(opt_prob.constraints_text_summary())
            opt_prob.resolve_constraints()
            # print(opt_prob.all_constraints_pass())
            opt_prob.optimize()
            # print(opt_prob.constraints_text_summary())
            # print(opt_prob.objectives_text_summary())
            if opt_prob.all_constraints_pass():
                modified_sequence = Seq(opt_prob.sequence)
            else:
                print("Warning: forbidden motifs could not be eliminated from mutatant %s." % mutant, file=sys.stderr)

        insert_sequence = modified_sequence[left_edge:right_edge]
        # print()
        # print(mutant)
        # print(insert_sequence)

        # TODO: account for the case where the regions on the edges of the mutated region are actually not mutated,
        #  for example where position 3 in the new codon is mutated, but positions 1 and 2 are not.
        #      or just don't worry about that...
        left_primer_overhang = ""
        right_primer_overhang = ""
        # INTERNAL_PRIMER_EXTENSION
        if len(insert_sequence) > 6:
            raise ValueError("Error, mutated region too large %s" % mutant)
        else:
            left_primer_template = str(seq[0:left_edge].reverse_complement())
            right_primer_template = str(seq[right_edge:])

            left_primer_overhang = str(
                (insert_sequence + seq[right_edge:right_edge + params.five_prime_extension]).reverse_complement())
            right_primer_overhang = str(seq[left_edge - params.five_prime_extension:left_edge] + insert_sequence)

        left_primer_dict = find_primer(left_primer_template, params.targetTm, params.minTm, params.maxTm,
                                       params.min_length, params.max_length, params.Tris, params.Mg, params.K,
                                       params.Na, params.dNTPs, params.primer_conc, params.gc_clamp)
        right_primer_dict = find_primer(right_primer_template, params.targetTm, params.minTm, params.maxTm,
                                        params.min_length, params.max_length, params.Tris, params.Mg, params.K,
                                        params.Na, params.dNTPs, params.primer_conc, params.gc_clamp)

        if left_primer_dict is None:
            print("could not find left primer for seq %s" % mutant, file=sys.stderr)
        elif right_primer_dict is None:
            print("could not find right primer for seq %s" % mutant, file=sys.stderr)
        else:
            left_primer = left_primer_overhang + left_primer_dict['seq']
            right_primer = right_primer_overhang + right_primer_dict['seq']

            new_primer_tms[left_primer] = left_primer_dict['tm']
            new_primer_tms[right_primer] = right_primer_dict['tm']
            new_primers.append((target_mutant_names[mutant] + "_F", left_primer))
            new_primers.append((target_mutant_names[mutant] + "_R", right_primer))

    # with open(args.output_assembly_file,'w') as out:
    #  for (construct,parts) in assembly.items():
    #    print("\t".join([construct] + parts), file = out)

    if params.output_format == "fasta":
        out_linestart = ">"
        out_separator = "\n"
    else:
        out_linestart = ""
        out_separator = "\t"

    for (name, primer) in new_primers:
        if (len(primer) < MIN_PRIMER_LENGTH) or (len(primer) > MAX_PRIMER_LENGTH):
            print("Warning: primer %s is outside of hard size bounds %d,%d" %
                  (str(primer), MIN_PRIMER_LENGTH, MAX_PRIMER_LENGTH), file=sys.stderr)
        print("%s%s%s%s" % (out_linestart, name, out_separator, primer), file=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", default=None,
                        help="Fasta of input sequence. Should contain only one nucleotide sequence.")
    parser.add_argument("-o", default=None, help="File for primer sequences")
    parser.add_argument("--target_mutations", default=None, required=True,
                        help="File indicating the desired mutations.")
    parser.add_argument("--output_assembly_file", default=None, required=True,
                        help="File to save the assembly specifications.")
    parser.add_argument("--target_species_taxid", default=316407, type=int,
                        help="Taxid of target species. For mutants, the most common codon encoding "
                             "the target amino acid will be selected.")
    parser.add_argument("--minTm", type=float, default=50,
                        help="Will report an error if a primer can't be found with a Tm of at least this.")
    parser.add_argument("--maxTm", type=float, default=80,
                        help="Will report an error if a primer can't be found with a Tm of at most this.")
    parser.add_argument("--targetTm", type=float, default=70, help="Will try to find a primer with a Tm clost to this.")
    parser.add_argument("--min_length", type=int, default=20,
                        help="minimum acceptable length for the core primer (without the overhangs attached)")
    parser.add_argument("--max_length", type=int, default=35,
                        help="maximum acceptable length for the core primer (without the overhangs attached)")
    parser.add_argument("--enzymes_to_avoid", type=str, default=None,
                        help="Names of restriction enzymes whose sites should be excluded "
                             "from any new mutants. Comma separated.")
    parser.add_argument("--motifs_to_avoid", type=str, default=None,
                        help="Subsequences that should be excluded from any new mutants. Comma separated.")
    parser.add_argument("--output_format", type=str, default="table", choices={"table", "fasta"},
                        help="format to write output file in")
    parser.add_argument('--gc_clamp', action="store_true", default=False)
    parser.add_argument('--five_prime_extension', type=int, default=7,
                        help="the number of homologous bases to include 5' to the mutated site")

    parser.add_argument("--Na", type=float, default=50, help="Sodium ion concentration [50]")
    parser.add_argument("--K", type=float, default=0, help="Potassium ion concentration [0]")
    parser.add_argument("--Tris", type=float, default=0, help="Tris ion concentration [0]")
    parser.add_argument("--Mg", type=float, default=100, help="Magnesium ion concentration [100]")
    parser.add_argument("--dNTPs", type=float, default=50, help="dNTP ion concentration [50]")
    parser.add_argument("--primer_conc", type=float, default=250, help="primer concentration [250]")

    # parser.add_argument("--overhang_5", type=str, default="", help="5' overhang to attach to all primers")
    # parser.add_argument("--overhang_3", type=str, default="", help="3' overhang to attach to all primers")

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

    if args.i is not None:
        input_handle.close()

    if args.o is not None:
        output_handle.close()
