"""
DOCUMENTATION AND CONFIGURATION -----------------------------------------------
"""
# STANDARD LIBRARIES
import sys
import re
from itertools import groupby
# THIRD PARTY PACKAGES
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils.MeltingTemp import make_table, DNA_NN4, DNA_NN2
import python_codon_tables as pct
from dnachisel import (
    DnaOptimizationProblem,
    AvoidPattern,
    EnzymeSitePattern,
    DnaNotationPattern,
    EnforceSequence,
    EnforceTranslation,
    CodonOptimize,
)
# APPLICATION MODULES
import parsers
import idt_primer_output, Generate_Nanohive_excel
# CONFIG
MIN_PRIMER_LENGTH = 15
MAX_PRIMER_LENGTH = 60
INTERNAL_PRIMER_EXTENSION = 10


# ALGORITHM -------------------------------------------------------------------
def design_end_to_end_primers(primer_order):
    """Pass"""
    return

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


def find_primer(seq, target_tm, min_tm, max_tm, min_len, max_len, Tris, Mg, K, Na, dNTPs, primer_conc, gc_clamp, correction):
    """
    finds a subsequence starting at the beginning of seq, with properties as specified by the rest of the parameters.
    if no sequence matching the constraints can be found, None is returned.
    """
    # mytable = make_table(oldtable=DNA_NN2, values={'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8)})
    best_primer = dict()
    for l in range(min_len, max_len + 1):
        sub = seq[0:l]
        tm = MeltingTemp.Tm_NN(sub, Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs, saltcorr=correction, dnac1=primer_conc/2,
                               dnac2=primer_conc/2, nn_table=DNA_NN4)
        # tm = primer3.calcTm(sub, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
        #                     max_nn_length=60, tm_method='santalucia', salt_corrections_method='santalucia')
        if (sub[-1] in "GC") or (not gc_clamp):
            if (tm >= min_tm) and (tm <= max_tm):
                # print('IF1', tm, sub)
                # if nothing in best_primer and temp is less than the target
                if ("seq" not in best_primer) and (target_tm+3 - tm >= 0):
                    best_primer['seq'] = sub
                    best_primer['tm'] = tm
                    # print('IF2', tm, sub)
                # elif abs(target_tm - tm) < abs(target_tm - best_primer['tm']):
                #     best_primer['seq'] = sub
                #     best_primer['tm'] = tm
                # elif less than the target temp and closer to the target temp compared to the previous best
                elif (target_tm+3 - tm >= 0) and (abs(target_tm - tm) < abs(target_tm - best_primer['tm'])):
                    best_primer['seq'] = sub
                    best_primer['tm'] = tm
                    # print('IF3', tm, sub)
    if len(best_primer) == 0:
        return None
    else:
        # print(best_primer['tm'], best_primer['seq'])
        return best_primer

def split_text(s):
    for k, g in groupby(s, str.isnumeric):
        yield ''.join(g)

def main(inpt, output, params):
    target_mutations_file = params.target_mutations
    seqs = list(SeqIO.parse(inpt, "fasta"))
    if len(seqs) != 1:
        raise ValueError("Error: Input fasta must contain only one sequence")
    seq_name = seqs[0].name
    seq = seqs[0].seq

    # if (len(seq) % 3) != 0:
    #     raise ValueError("Error: Input sequence must have length divisible by 3")

    # mutations_file_lines = [x for x in parsers.SimpleParser(target_mutations_file)]
    # target_mutant_names = dict()
    # target_mutants = dict()
    # for line in mutations_file_lines:
    #     parts = [x.strip() for x in line.split("\t")]

    #     if len(parts) == 1:
    #         target_mutants[parts[0]] = parse_target_mutations(parts[0])
    #         target_mutant_names[parts[0]] = seq_name + "_" + parts[0]
    #     elif len(parts) == 2:
    #         target_mutants[parts[1]] = parse_target_mutations(parts[1])
    #         target_mutant_names[parts[1]] = parts[0]
    #     else:
    #         print("Warning, too many fields in design file line: %s" % line, file=sys.stderr)
    #

    if str(target_mutations_file).endswith('.xlsx'):
        df = pd.read_excel(target_mutations_file, sheet_name='Mutations')
        target_mutant_names = dict()
        target_mutants = dict()
        count = int(df.iloc[0, 1])  # Start with the number in the input index ID example
        primer_type = []
        list_of_complex_mutant_codes = [x for x in list(df.iloc[2:, 0]) if str(x) != 'nan']
        for code in list_of_complex_mutant_codes:
            code_list = list(split_text(code))
            base_code = code_list[0] + code_list[1]
            targets = ''.join(code_list[2:])
            idx = 0
            while idx < len(targets):
                character = targets[idx]
                if character not in ['(', ')', '#', '!']:  # run normally
                    mutation_code = base_code + character
                    target_mutants[mutation_code] = parse_target_mutations(mutation_code)
                    target_mutant_names[mutation_code] = str(df.columns[1])+'-'+str(count).zfill(4)
                    primer_type.append('Single')
                    count += 1; idx += 1
                elif character == '(':  # pull the next three from targets to get nt mutation code
                    mutation_code = base_code + targets[idx+1 : idx+4]
                    target_mutants[mutation_code] = parse_target_mutations(mutation_code)
                    target_mutant_names[mutation_code] = str(df.columns[1])+ '-' + str(count).zfill(4)
                    primer_type.append('Degenerate')
                    count += 1; idx += 5
                elif character == '!':  # run for smart primers; Tang L, et al. Construction of “small-intelligent” ... Biotechniques. 2012
                    for nt_code, id_letter in zip(['NDT','VMA','ATG','TGG'], ['a','b','c','d']):
                        mutation_code = base_code + nt_code
                        target_mutants[mutation_code+'smart'] = parse_target_mutations(mutation_code)
                        target_mutant_names[mutation_code+'smart'] = str(df.columns[1])+'-'+str(count).zfill(4)+id_letter
                        primer_type.append('SMART')
                    count += 1; idx += 1
                elif character == '#':  # run for trick primers; Kille, S. et al. Reducing codon redundancy... ACS Synth. Biol. (2013).
                    for nt_code, id_letter in zip(['NDT','VHG','TGG'], ['a','b','c']):
                        mutation_code = base_code + nt_code
                        target_mutants[mutation_code+'trick'] = parse_target_mutations(mutation_code)
                        target_mutant_names[mutation_code+'trick'] = str(df.columns[1])+'-'+str(count).zfill(4)+id_letter
                        primer_type.append('TRICK')
                    count += 1; idx += 1
                else:
                    print('error with mutation codes', code)
                    idx += 1
    else:
        print('Please provide either an excel file (.xlsx)')

    # target_mutants = {x: parse_target_mutations(x) for x in parsers.simple_parser(target_mutations_file)}
    # #dict where key is mutation specification, and value is the parsed mutation specification
    # (a lists [mutation number] of dicts {sourceaa, targetaa, position} )

    # new_primers = dict() #key is sequence, value is list of target_mutations that the sequence corresponds to.
    new_primer_names = dict()
    assembly = dict()
    new_primer_tms = dict()
    products = dict()

    # returns a dict where keys are AA and values are dicts where keys are codons and values are normalized frequency.
    aa_to_codons = pct.get_codons_table(
        params.target_species_taxid)

    if len(aa_to_codons) == 0:
        raise ValueError("Error: cannot find codon table corresponding to taxid:" % params.target_species_taxid)
    codon_to_aa = invert_codon_table(aa_to_codons)

    sorted_codon_table = sort_codon_table(aa_to_codons)

    counters = {'left': 0, 'right': 0}

    # mutation_space=list() #for codon optimization

    for (mutant, mutations), primer_type_temp in zip(target_mutants.items(), primer_type):
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
            # mutation_space.append( MutationSpace.MutationChoice( (nuc_pos,nuc_pos+3), \
            # { x[0] for x in sorted_codon_table[mutations[i]['targetaa']] }) )
            modified_sequence = modified_sequence[0:nuc_pos] + new_codon + modified_sequence[nuc_pos + 3:]
        modified_sequence = Seq(modified_sequence)

        # TODO: right here is where we'd check modified sequence to make sure it
        #  doesn't have any forbidden restriction sites.
        if (params.enzymes_to_avoid is not None) or (params.motifs_to_avoid is not None):
            constraints = list()
            constraints.append(EnforceSequence(sequence=modified_sequence[0:left_edge], location=(0, left_edge, 1)))
            constraints.append(EnforceSequence(sequence=modified_sequence[right_edge:],
                                               location=(right_edge, len(modified_sequence), 1)))
            constraints.append(EnforceTranslation(genetic_table='Standard', location=(0, len(modified_sequence), 1),
                                                  start_codon='keep'))
            # TODO: Handle if genetic table is not Standard, also handle if this is a subsequence
            #  of a larger sequence, possibly handled via annotations.
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
        if len(insert_sequence) <= 6:
            left_primer_overhang = str(insert_sequence.reverse_complement())
            #right_primer_overhang = str(insert_sequence)                                                                # Added by TJM
        else:
            mid = len(insert_sequence) // 2
            left_primer_overhang = str(insert_sequence[0:mid].reverse_complement())
            right_primer_overhang = str(insert_sequence[mid:])

        left_primer_template = str(seq[0:left_edge].reverse_complement())

        if params.internal_mutation and (len(insert_sequence) <= 6):
            right_primer_template = str(seq[right_edge + INTERNAL_PRIMER_EXTENSION:])
        else:
            right_primer_template = str(seq[right_edge:])

        # TODO:account for INTERNAL_PRIMER_EXTENSION in bounds for left_primer
        left_ballast = len(left_primer_overhang)
        if params.internal_mutation and (len(insert_sequence) <= 6):
            left_ballast += INTERNAL_PRIMER_EXTENSION
        left_primer_dict = find_primer(left_primer_template, params.targetTm, params.minTm, params.maxTm,
                                       max(10, params.min_length - left_ballast),
                                       max(10, params.max_length - left_ballast),
                                       params.Tris, params.Mg, params.K, params.Na, params.dNTPs, params.primer_conc,
                                       params.gc_clamp, params.correction)
        if left_primer_dict is not None:
            print('left ', left_primer_dict['tm'], left_primer_dict['seq'])
        # print(right_primer_overhang, max(10, params.min_length - len(right_primer_overhang)))
        right_primer_dict = find_primer(right_primer_template, params.targetTm, params.minTm, params.maxTm,
                                        # max(10, params.min_length - len(right_primer_overhang)),
                                        # max(10, params.max_length - len(right_primer_overhang)), params.Tris, params.Mg,
                                        max(10, params.min_length - 3),
                                        max(10, params.max_length - 3), params.Tris, params.Mg,
                                        params.K, params.Na, params.dNTPs, params.primer_conc,
                                        params.gc_clamp, params.correction)
        if right_primer_dict is not None:
            print('right', right_primer_dict['tm'], right_primer_dict['seq'])

        if left_primer_dict is None:
            print("could not find left primer for seq %s" % mutant, file=sys.stderr)
        elif right_primer_dict is None:
            print("could not find right primer for seq %s" % mutant, file=sys.stderr)
        else:
            if params.internal_mutation and (len(insert_sequence) <= 6):
                left_primer = str(seq[
                                  right_edge:right_edge + INTERNAL_PRIMER_EXTENSION].reverse_complement()) + \
                              left_primer_overhang + \
                              left_primer_dict['seq']
            else:
                left_primer = left_primer_overhang + left_primer_dict['seq']

            right_primer = right_primer_overhang + right_primer_dict['seq']

            new_primer_tms[left_primer] = left_primer_dict['tm']
            new_primer_tms[right_primer] = right_primer_dict['tm']
            # print(left_primer)
            # print(right_primer)

            if left_primer not in new_primer_names:
                #     new_primers[left_primer] = list()
                counters['left'] += 1
                new_primer_names[left_primer] = seq_name + "_" + "tm-" + (
                        "%.1f" % new_primer_tms[left_primer]) + "_" + "l_" + str(counters['left'])
            if right_primer not in new_primer_names:
                #      new_primers[right_primer] = list()
                counters['right'] += 1
                new_primer_names[right_primer] = seq_name + "_" + "tm-" + (
                        "%.1f" % new_primer_tms[right_primer]) + "_" + "r_" + str(counters['right'])

            # new_primers[right_primer].append(mutant)
            # new_primers[left_primer].append(mutant)

            construct_name = target_mutant_names[mutant]
            assembly[construct_name] = [seq_name, new_primer_names[left_primer], new_primer_names[right_primer], mutant, left_primer_dict['tm'], params.targetTm, primer_type_temp]
            products[construct_name] = modified_sequence


    # with open(params.output_assembly_file, 'w') as out:
    #     for (construct, parts) in assembly.items():
    #         print("\t".join([construct] + parts), file=out)

    if params.output_format == "fasta":
        out_linestart = ">"
        out_separator = "\n"
    else:
        out_linestart = ""
        out_separator = "\t"

    new_primers = list()
    for primer in new_primer_names:
        print(primer)
        if (len(primer) < MIN_PRIMER_LENGTH) or (len(primer) > MAX_PRIMER_LENGTH):
            print("Warning: primer %s is outside of hard size bounds %d,%d" % (
                str(primer), MIN_PRIMER_LENGTH, MAX_PRIMER_LENGTH), file=sys.stderr)
        print("%s%s%s%s" % (out_linestart, new_primer_names[primer], out_separator, primer), file=output)
        new_primers.append((new_primer_names[primer], Seq(primer)))

    if params.output_products_file is not None:
        with open(params.output_products_file, "w") as outfile:
            for n, s in products.items():
                print(">" + n + "\n" + s, file=outfile)


    # write out primers to file formatted for IDT upload
    # read_primers = idt_primer_output.read_input_file(params.o)
    idt_primer_output.output_idt(params.o, new_primers, '')
    Generate_Nanohive_excel.output_nanohive_excel(params.paths, assembly, params.final_rxn_volume,
                                                  params.final_primer_conc, params.final_template_conc)
    Generate_Nanohive_excel.output_nanohive_excel2(params.paths, assembly, params.final_rxn_volume,
                                                   params.final_primer_conc, params.final_template_conc)
    Generate_Nanohive_excel.output_biomek_excel_files(params.paths, assembly, params.final_rxn_volume,
                                                      params.final_primer_conc)
    Generate_Nanohive_excel.output_biomek_excel_file_all_primers(params.paths, assembly, params.final_rxn_volume,
                                                      params.final_primer_conc)
