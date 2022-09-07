"""Contains the Primer Object

Description.

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
import subprocess
import collections
# Third Party Packages
import regex
from Bio.SeqUtils import MeltingTemp
# Local Modules
from .dna import DNA
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASS -----------------------------------------------------------------------
class Primer(DNA):
    """Representation of a primer"""
    def __init__(self, sequence, name='generic_primer', project='Generic',
        purpose='Unspecified', method='Unspecified'):
        """Primer constructor function"""
        super().__init__(sequence=sequence, name=name)
        self.purpose = purpose
        self.method = method
        self.project = project
        return

    @property
    def tm(self):
        "Tm of the primer if it were to anneal to its reverse complement"
        return DNA.calculate_nn_tm(self.sequence)

    def anneal(self, template, required_three_prime_match=10, required_tm=40,
        additional_five_prime_match=15, allow_isolated_mismatch=True):
        """
        Returns the binding index, mismatch index, and direction of the
        primer on the specified template.
        """
        if template.length > 1_000_000:
            required_three_prime_match = 15
            required_tm = 40
            allow_isolated_mismatach=False
        primer_binding_sites = find_binding_sites(template)
        if primer_binding_sites:
            annealed_primers = determine_primer_match_index(template, primer_binding_sites)
        else:
            raise Exception('No binding sites could be found on the inputted template.')

        def find_binding_sites(self, template):
            """Return list of 3' indices for all identified primer binding sites"""
            forward_search = regex.finditer(
                self.sequence[-required_three_prime_match:]+r'{e<=1}',
                template.sequence,
                overlapped=True)
            reverse_search = regex.finditer(
                self.sequence[-required_three_prime_match:]+r'{e<=1}',
                template.reverse_complement,
                overlapped=True)
            binding_sites = [hit.end()-1 for hit in forward_search] + [-hit.end() for hit in reverse_search]
            return binding_sites

        def determine_primer_match_index(self, template, binding_sites):
            """
            Returns list of boolean values where each value represents if the
            nucleotide at the respective primer index anneals to the template
            """
            for binding_site in binding_sites:
                if binding_site >= 0:
                    site = binding_site
                    template_sequence = template.sequence
                else:
                    site = -(binding_site+1)
                    template_sequence = template.reverse_complement
                for nucleotide in reversed(self.sequence):
                    if nucleotide == template_sequence[site]:
                        site -= 1
                    else:
                        mismatch_index = (binding_site-self.length, site)
                        yield AnnealedPrimer(self, template, binding_site, mismatch_index)
                yield AnnealedPrimer(self, template, binding_site, None)
        
        return list(annealed_primers)

    def blast(self, template):
        """Blasts the primer sequence against a template sequence to search for
        binding sites."""
        TEMP_DIR = ROOT_DIR / 'static' / 'temp'
        template_fasta = TEMP_DIR / 'template.fasta'
        template_fasta.touch()
        with template_fasta.open(mode='w') as template_file:
            template_file.write(f'>{template.name}\n{template.sequence}')
        primer_fasta = TEMP_DIR / 'primer.fasta'
        primer_fasta.touch()
        with primer_fasta.open(mode='w') as primer_file:
            primer_file.write(f'>{self.name}\n{self.sequence}')
        template_db = TEMP_DIR / 'template_db'
        subprocess.run(f'makeblastdb -dbtype nucl -in {template_fasta}  -title {template.name} -out {template_db}', shell=True, check=True)
        hits_csv = TEMP_DIR / 'hits.csv'
        subprocess.run(f'blastn -db {template_db} -query {primer_fasta} -evalue 1000 -outfmt "10 qseqid qlen sseqid slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand" -out {hits_csv}', shell=True, check=True)
        subprocess.run(f'rm -rf {TEMP_DIR}/*db*', shell=True, check=True)
        with hits_csv.open(mode='r') as file:
            string_output = 'query_sequence_id,query_sequence_length,subject_sequence_id,subject_sequence_length,query_alignment_start,query_alignment_end,subject_alignment_start,subject_alignment_end,aligned_query_sequence,aligned_subject_sequence,e-value,bitscore,raw_score,alignment_length,percent_identical_matches,identical_matches,mismatches,positive-scoring_matches,gap_openings,gaps,percentage_positive-scoring_matches,subject_strand\n'+file.read()
        template_fasta.unlink()
        primer_fasta.unlink()
        hits_csv.unlink()
        return string_output

    @classmethod
    def design_primer(cls):
        """Design a single primer based on the input parameters"""
        return

    @classmethod
    def design_primer_pair(cls):
        """Design a primer_pair based on the input parameters"""
        return


class AnnealedPrimer(Primer):
    """Representation of one nucleotide bound to another"""
    def __init__(self, primer, template, binding_site, mismatch_index):
        self.primer = primer
        self.template = template
        # integer representing the indec of the 3' nucleotide
        self.binding_site = binding_site

        self.mismatch_index = mismatch_index
        return

    @property
    def tm1(self):
        """Returns the tm of the primer on the bound template"""
        if self.binding_site >=0:
            tm = DNA.calculate_nn_tm(
                seq=self.primer.sequence,
                c_seq=self.template.complement[self.binding_site-self.length-1:self.binding_site+1],
                shift=1)
        else:
            tm = DNA.calculate_nn_tm(
                seq=self.primer.sequence,
                c_seq=self.sequence[self.binding_site+self.length-1:self.binding_site+1:-1],
            )
        return tm

    @property
    def tm2(self):
        """Returns the tm of the primer on the bound template"""
        if self.binding_site >=0:
            tm = DNA.calculate_nn_tm(
                seq=self.primer.sequence,
                c_seq=self.primer.complement + self.template.complement[self.binding_site+1])
        else:
            tm = DNA.calculate_nn_tm(
                seq=self.primer.sequence,
                c_seq=self.primer.complement + self.template.reverse()[self.binding_site-1])
        return tm

    @staticmethod
    def design_overlapping_primers(primer_order, codon_index, codon): # incomplete
        """
        Generates and selects the best primer pair for the input mutations

        Args:
            primer_order
            codon_index
            codon

        Returns:
            primer_pair
        """
        def main():
            for variant_code in primer_order.variant_codes:
                forward_starting_primer = generate_forward_starter_primer(codon_index, codon)
                forward_primers = generate_primer_candidates(forward_starting_primer)
                reverse_starting_primer = generate_reverse_starter_primer(codon_index, codon)
                reverse_primers = generate_primer_candidates(reverse_starting_primer)
                forward_primer, reverse_primer = select_best_primer_pair(forward_primers, reverse_primers)
                primer_pair = (forward_primer, reverse_primer)
                yield primer_pair

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

        return main()

class PrimerPair:
    """Placeholder"""

    def __init__(self, annealed_primer_1, annealed_primer_2):
        self.f_primer = annealed_primer_1
        self.r_primer = annealed_primer_2
