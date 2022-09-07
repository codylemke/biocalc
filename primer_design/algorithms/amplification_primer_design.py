"""

"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'


# ALGORITHM -------------------------------------------------------------------
def design_amplification_primers(primer_order):
    """
    Constructs minimal viable primers for mutagenesis using the defined primer parameters.
    """
    main(primer_order)

    def main(primer_order):
        """Algorithm logic"""
        for variant_code in primer_order.variant_codes:
            forward_starting_primer = generate_minimal_starter_primer(codon_index, codon)
            forward_primers = generate_primer_candidates(forward_starting_primer)
            reverse_starting_primer = generate_reverse_starter_primer(codon_index, codon)
            reverse_primers = generate_primer_candidates(reverse_starting_primer)
            forward_primer, reverse_primer = select_best_primer_pair(forward_primers, reverse_primers)
            primer_pair = (forward_primer, reverse_primer)