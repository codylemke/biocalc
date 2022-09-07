"""

-------------------------------------------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

# STANDARD LIBRARIES
import re
# APPLICATION MODULES
from primer_design.objects import Primer


# ALGORITHM -------------------------------------------------------------------
def design_overlapping_primers(primer_order, codon_index, codon):
    """
    Constructs minimal viable primers for mutagenesis using the defined primer parameters.
    """
    for variant_code in primer_order.variant_codes:
        forward_starting_primer = generate_forward_starter_primer(codon_index, codon)
        forward_primers = generate_primer_candidates(forward_starting_primer)
        reverse_starting_primer = generate_reverse_starter_primer(codon_index, codon)
        reverse_primers = generate_primer_candidates(reverse_starting_primer)
        forward_primer, reverse_primer = select_best_primer_pair(forward_primers, reverse_primers)
        primer_pair = (forward_primer, reverse_primer)

    def generate_forward_starter_primer(codon_index, codon):
        """Returns a starter forward primer using the minimal edge length constraints"""
        primer_name = f'{primer_order.template.name}_{variant_code}_F'
        primer_index = [codon_index[0] - primer_order.MIN_FIVE_PRIME_EDGE_LENGTH, codon_index[1] + primer_order.MIN_THREE_PRIME_EDGE_LENGTH]
        forward_starter_primer = Primer(
            name=primer_name,
            sequence=f'{primer_order.template.sequence[primer_index[0]:codon_index[0]]}{primer_order.template.sequence[codon_index[0]:codon_index[1]]}{primer_order.template.sequence[codon_index[1]:primer_index[1]]}',
            project=primer_order.project,
            purpose=primer_order.purpose,
            method=primer_order.method,
            index=primer_index,
            mismatch_index=codon_index,
            template=primer_order.template.sequence,
            order=primer_order
        )
        return forward_starter_primer

    def generate_reverse_starter_primer(codon_index, codon):
        """Returns a starter reverse primer using the minimal edge length constraints"""
        primer_name = f'{primer_order.template.name}_{variant_code}_R'
        rc_codon_index = (primer_order.template.length - codon_index[1] - 1, primer_order.template.length - codon_index[1] - 1)
        primer_index = [rc_codon_index[0] - primer_order.MIN_FIVE_PRIME_EDGE_LENGTH, rc_codon_index[1] + primer_order.MIN_THREE_PRIME_EDGE_LENGTH]
        reverse_starter_primer = Primer(
            name=primer_name,
            sequence=f'{primer_order.template.complement()[primer_index[0]:rc_codon_index[0]]}{primer_order.template.complement()[rc_codon_index[0]:rc_codon_index[1]]}{primer_order.template.complement()[rc_codon_index[1]:primer_index[1]]}',
            project=primer_order.project,
            purpose=primer_order.purpose,
            method=primer_order.method,
            index=primer_index,
            mismatch_index=rc_codon_index,
            template=primer_order.template.complement(),
            order=primer_order
        )
        return reverse_starter_primer

    def generate_primer_candidates(primer):
        """Returns list of all valid primers given the constraints"""
        # Generate list of all primers within the length constrains
        primers = list()
        while primer.five_prime_edge_length < primer_order.MAX_FIVE_PRIME_EDGE_LENGTH and primer.length < primer_order.MAX_PRIMER_LENGTH:
            primer.index[0] -= 1
            primer.save()
            if primer_order.validate_primer(primer):
                primers.append(primer)
        for primer in primers:
            while primer.three_prime_edge_length < primer_order.MAX_THREE_PRIME_EDGE_LENGTH and primer.length < primer_order.MAX_PRIMER_LENGTH:
                primer.index[1] += 1
                primer.save()
                if primer_order.validate_primer(primer):
                    primers.append(primer)
        # Remove duplicate primers
        primer_filter = list()
        unique_primers = list()
        for primer in primers:
            if primer.sequence not in primer_filter:
                primer_filter.append(primer.sequence)
                unique_primers.append(primer)
        # Score unique primers
        primer_candidates = list()
        for primer in unique_primers:
            scored_primer = score_primer(primer)
            primer_candidates.append(scored_primer)

        def score_primer(primer_order, primer):
            """Returns the input primer with primer scores based on optimal parameters"""
            five_prime_edge_length_range = abs(primer_order.MAX_FIVE_PRIME_EDGE_LENGTH - primer_order.MIN_FIVE_PRIME_EDGE_LENGTH)
            three_prime_edge_length_range = abs(primer_order.MAX_THREE_PRIME_EDGE_LENGTH - primer_order.MIN_THREE_PRIME_EDGE_LENGTH)
            primer_length_range = abs(primer_order.MAX_PRIMER_LENGTH - primer_order.MIN_PRIMER_LENGTH)
            five_prime_edge_tm_range = abs(primer_order.MAX_FIVE_PRIME_EDGE_TM - primer_order.MIN_FIVE_PRIME_EDGE_TM)
            three_prime_edge_tm_range = abs(primer_order.MAX_THREE_PRIME_EDGE_TM - primer_order.MIN_THREE_PRIME_EDGE_TM)
            primer_tm1_range = abs(primer_order.MAX_PRIMER_TM1 - primer_order.MIN_PRIMER_TM1)
            primer_tm2_range = abs(primer_order.MAX_PRIMER_TM2 - primer_order.MIN_PRIMER_TM2)
            primer_gc_range = abs(primer_order.MAX_PRIMER_GC - primer_order.MIN_PRIMER_GC)
            primer.five_prime_edge_length_score = (
                abs(primer_order.OPTIMAL_FIVE_PRIME_EDGE_LENGTH - primer.five_prime_edge_length)
                / five_prime_edge_length_range
                * primer_order.FIVE_PRIME_EDGE_LENGTH_BIAS
            )
            primer.three_prime_edge_length_score = (
                abs(primer_order.OPTIMAL_THREE_PRIME_EDGE_LENGTH - primer.three_prime_edge_length)
                / three_prime_edge_length_range
                * primer_order.THREE_PRIME_EDGE_LENGTH_BIAS
            )
            primer.length_score = (
                abs(primer_order.OPTIMAL_PRIMER_LENGTH - primer.length)
                / primer_length_range
                * primer_order.PRIMER_LENGTH_BIAS
            )
            primer.five_prime_edge_tm_score = (
                abs(primer_order.OPTIMAL_FIVE_PRIME_EDGE_TM - primer.five_prime_edge_tm)
                / five_prime_edge_tm_range
                * primer_order.FIVE_PRIME_EDGE_TM_BIAS
            )
            primer.three_prime_edge_tm_score = (
                abs(primer_order.OPTIMAL_THREE_PRIME_EDGE_TM - primer.three_prime_edge_tm)
                / three_prime_edge_tm_range
                * primer_order.THREE_PRIME_EDGE_TM_BIAS
            )
            primer.tm_score = (
                abs(primer_order.OPTIMAL_PRIMER_TM1 - primer.tm)
                / primer_tm1_range
                * primer_order.PRIMER_TM1_BIAS
            )
            primer.tm2_score = (
                abs(primer_order.OPTIMAL_PRIMER_TM2 - primer.tm2)
                / primer_tm2_range
                * primer_order.PRIMER_TM2_BIAS
            )
            primer.gc_score = (
                abs(primer_order.OPTIMAL_PRIMER_GC - primer.gc_content)
                / primer_gc_range
                * primer_order.PRIMER_GC_BIAS
            )
            primer.score = (
                primer.five_prime_edge_length_score
                + primer.three_prime_edge_length_score
                + primer.length_score
                + primer.five_prime_edge_tm_score
                + primer.three_prime_edge_tm_score
                + primer.tm1_score
                + primer.tm2_score
                + primer.gc_score
            )
            return scored_primer
        
        return primer_candidates          

    def select_best_primer_pair(forward_primers, reverse_primers):
        """Returns the best primer pair based on scoring metrics"""
        # Unfinished
        primer_pairs = list()
        for forward_primer in forward_primers:
            for reverse_primer in reverse_primers:
                primer_pair = (forward_primer, reverse_primer)
        best_primer_pair = None
        for primer_pair in primer_pairs:
            if primer_pair.score > best_primer_pair.score:
                best_primer_pair = primer_pair
        return best_primer_pair

    return primer_pair

