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
    def tmi(self):
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
    def tmf(self):
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