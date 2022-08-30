"""Contains the DNA object

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
import re
import collections
# Third Party Packages
import numpy as np
from Bio.SeqUtils import MeltingTemp
# Local Modules
from .nucleotide import Nucleotide
from .protein import Protein
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class DNA(Nucleotide):
    """Representation of DNA"""
    
    def __init__(
        self, sequence, name='generic_dna', dna_type='generic'):
        """Constructor Function"""
        super().__init__(sequence=sequence, name=name, nucleotide_type='dna')
        self.dna_type = dna_type
        return

    def annotate(self):
        """
        Loops through a database of DNA features saerching the sequence to
        see if the feature is present. Saves a list of features and respective
        indeces as an attribute of the DNA object
        """
        return

    def find_cut_sites(self, enzymes):
        """
        Scans the sequence looking for restriction enzyme binding sites and
        returns the indices of the binding sites and their respective cut sites
        """
        return # [[index, recognition_site, cut_site]

    def find_orfs(self):
        """
        Scans the sequence looking for start codons then scans for the first
        codon in frame with the start cofon and returns the orf indices
        """
        Orf = collections.namedtuple('Orf', ['strand', 'start', 'stop', 'sequence'])
        orfs = list()
        for start_codon in re.finditer('ATG', self.sequence):
            for index in range(start_codon.end(), self.length, 3):
                codon = self.sequence[index:index+3]
                if codon in ['TAA', 'TGA', 'TAG']:
                    orf = Orf('+', start_codon.start(), index+3, self.sequence[start_codon.start():index+3])
                    orfs.append(orf)
                    break
        for start_codon in re.finditer('ATG', self.reverse_complement):
            for index in range(start_codon.end(), self.length, 3):
                codon = self.sequence[index:index+3]
                if codon in ['TAA', 'TGA', 'TAG']:
                    orf = Orf('-', start_codon.start(), index+3, self.reverse_complement[start_codon.start():index+3])
                    orfs.append(orf)
                    break
        return orfs

    @staticmethod
    def calculate_nn_tm(
        seq, check=True, strict=False, c_seq=None, shift=0, dnac1=25, dnac2=25,
        selfcomp=False, Na=50, K=0, Tris=0, Mg=1.5, dNTPs=0.6, saltcorr=7):
        """Calculates melting temperature using BioPython's most up to date
        nearest-neighbor methods."""
        tm = int(MeltingTemp.Tm_NN(
            seq=seq,
            check=check,
            strict=strict, # Note: thermodynamic data is not available for all mismatches which requires the strict parameter to be set to False to prevent failing in these scenarios
            c_seq=c_seq,
            shift=shift,
            nn_table=MeltingTemp.DNA_NN4, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            tmm_table=MeltingTemp.DNA_TMM1, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            imm_table=MeltingTemp.DNA_IMM1, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            de_table=MeltingTemp.DNA_DE1, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            dnac1=dnac1, # primer concentration
            dnac2=dnac2, # template_concentration
            selfcomp=selfcomp,
            Na=Na,
            K=K,
            Tris=Tris,
            Mg=Mg,
            dNTPs=dNTPs,
            saltcorr=saltcorr, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        ))
        return tm

    @staticmethod
    def calculate_wallace_tm(seq, check=True, strict=False):
        """
        Calculates melting temperature using BioPython's wallace method.
        """
        tm = int(MeltingTemp.Tm_Wallace(
            seq=seq,
            check=check,
            strict=strict, # Note: thermodynamic data is not available for all mismatches which requires the strict parameter to be set to False to prevent failing in these scenarios
        ))
        return tm

    @staticmethod
    def calculate_gc_tm(
        seq, check=True, strict=False, valueset=7, userset=None, Na=50, K=0,
        Tris=0, Mg=1.5, dNTPs=0.6, saltcorr=7, mismatch=True):
        """
        Calculates melting temperature using BioPython's gc method.
        """
        tm = int(MeltingTemp.Tm_GC(
            seq=seq,
            check=check,
            strict=strict, # Note: thermodynamic data is not available for all mismatches which requires the strict parameter to be set to False to prevent failing in these scenarios
            valueset=valueset,
            userset=userset,
            Na=Na,
            K=K,
            Tris=Tris,
            Mg=Mg,
            dNTPs=dNTPs,
            saltcorr=saltcorr, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            mismatch=mismatch
        ))
        return tm

    def anneal_oligo(self, primer):
        """Finds the template index based on the sequence"""
        # Find if the last 13 bases of the primer anneal to the template
        three_prime_seed_length = 13
        try:
            binding_sites = [index.start() for index in re.finditer(primer.sequence[-three_prime_seed_length:], self.sequence)]
        except ValueError as err:
            err('Annealing site not found on template sequence.')
        else:
            for binding_site in binding_sites:
                binding_index = [binding_site, binding_site + three_prime_seed_length]
                # Find how much of the 3' End of the primer anneals to the template
                while True:
                    three_prime_seed_length += 1
                    try:
                        binding_index[0] = self.sequence.index(primer.sequence[-three_prime_seed_length:])
                    except ValueError:
                        three_prime_seed_length -= 1
                        annealed_three_prime_sequence = self.sequence[-three_prime_seed_length:]
                        primer_mismatch_index = [0, primer.sequence.index(annealed_three_prime_sequence)]
                        break
                    else:
                        if len(self.sequence[-three_prime_seed_length:]) == len(primer.sequence):
                            primer_mismatch_index = None
                            return binding_index, primer_mismatch_index
                # Find if the 5' end of the primer binds to the template
                five_prime_seed_length = 13
                try:
                    self.sequence.index(primer.sequence[:five_prime_seed_length])
                except ValueError:
                    return binding_index, primer_mismatch_index
                else:
                # Find how much of the 5' end of the primer anneals to the template
                    while True:
                        five_prime_seed_length += 1
                        try:
                            self.sequence.index(primer.sequence[:five_prime_seed_length])
                        except ValueError:
                            five_prime_seed_length -= 1
                            annealed_five_prime_sequence = primer.sequence[:five_prime_seed_length]
                            binding_index[0] = self.sequence.index(primer.sequence[:five_prime_seed_length])
                            primer_mismatch_index[0] = binding_index[0] + len(annealed_five_prime_sequence)
                            return binding_index, primer_mismatch_index

    def pcr_amplify(self, primer_1=0, primer_2=None, project='generic'):
        """Finds where the primers bind in the sequence and returns the
        sequence that would result in PCR amplification"""
        from .fragment import Fragment
        PrimerPair = collections.namedtuple('PrimerPair', ['f_primer', 'r_primer'])
        # Polish the primer.blast() (anneal) function
        # Create AnnealedPrimer object based on blast output
        binding_site_1 = primer_1.blast(self)
        binding_site_2 = primer_2.blast(self)
        if binding_site_1.strand == '+' and binding_site_2.strand == '-':
            primer_pair = PrimerPair(binding_site_1, binding_site_2)
        elif binding_site_1.strand == '-' and binding_site_2.strand == '+':
            primer_pair = PrimerPair(binding_site_2, binding_site_1)
        fragment = Fragment(sequence=primer_1.sequence+self.sequence[binding_site_1])
        # In Development
        return # Fragement

    def design_amplification_primers(self, index_1=0, index_2=None, project='generic', left_extension=None, right_extension=None):
        """Returns a pair of primer objects dependent on the input indices for
        amplification"""
        Constraints = collections.namedtuple('Constraints', [])
        return # (f_primner, r_primer)

    def find_golden_gate_overhangs(self, overhangs):
        """Returns tuples of the index an overhang was found and the overhang
        that was found at that index"""
        return # [index, overhang, direction]

    def transcribe(self, strand='+', index_1=None, index_2=None):
        """Returns an RNA object with a sequence transcribed from the given
        sequence"""
        if strand == '+':
            dna = np.array(list(self.sequence[index_1:index_2]))
        elif strand == '-':
            dna = np.array(reversed(self.sequence[index_1:index_2]))
        rna = np.where(dna == 'T', 'U', dna)
        return ''.join(rna)

    def generate_gge_fragment(self, organism, module):
        """Returns a Fragment object with a sequence that allows it to be used
        in the Golden Gate Express platform"""
        from .fragment import Fragment

        def main():
            modules = parse_module(module)
            prepared_sequence = prepare_sequence(modules)
            sequence = append_gge_adapters(modules, prepared_sequence)
            fragment = Fragment(name=self.name, sequence=sequence)
            return fragment

        def parse_module(module):
            """Returns a named tuple representing the desired modules that the
            sequence is to span."""
            Modules = collections.namedtuple('Module', 'left right')
            if len(module_list := list(str(module))) == 1:
                modules = Modules(int(module_list[0]), int(module_list[0]))
            elif len(module_list) == 2:
                modules = Modules(int(module_list[0]), int(module_list[1]))
            else:
                raise ValueError('An invalid module number was selected')
            return modules

        def prepare_sequence(modules):
            """Alters the sequence if necessary to prepare it for adding
            adapters"""
            sequence = self.sequence
            tis = {
                    'e_coli': 'AGGAGAGCAGCTATG',
                    'e_coli_enhanced': 'AGGAGAGCAGCTATGCAGCTT',
                    'yeast': 'AGGAAAAAAATGTCT',
                    'mammalian': 'AGGAGCCACCATGGGC',
                    'plant': 'ACAACAATGGCT',
                    'microalga': 'GCCAAGATGGCG',
                }
            if modules.right == 3:
                if organism == 'e_coli':
                    rbs_agga_index = sequence.rindex('AGGA', start=-60)
                    sequence = sequence[:rbs_agga_index]
            if modules.left == 4:
                if sequence[:3] == 'ATG':
                    sequence = sequence[3:]
                try:
                    sequence = tis[organism]+sequence
                except KeyError as err:
                    raise KeyError('The organism entered is invalid.') from err
            if modules.right == 4:
                if sequence[-3:] in ['TAA', 'TGA', 'TAG']:
                    sequence = sequence[:-3]
                sequence = sequence+'AGTAGTG'
            if modules.left == 5:
                if sequence[:3] == 'ATG':
                    sequence = sequence[3:]
                sequence = 'AGTGGT'+sequence
            if modules.right == 5:
                if sequence[-3:] in ['TAA', 'TGA', 'TAG']:
                    sequence = sequence[:-3]
                sequence = sequence+'GGTAGCAGC'
            if modules.left == 6:
                if sequence[:3] == 'ATG':
                    sequence = sequence[3:]
            if modules.right == 6:
                if sequence[-3:] != 'TAA':
                    if sequence[-3:] == 'TGA' or sequence[-3:] == 'TAG':
                        sequence = sequence[:-3]+'TAA'
                    else:
                        sequence = sequence+'TAA'
            return sequence

        def append_gge_adapters(modules, sequence):
            """Returns the sequence of the portion of the overhang that should be
            appended to the sequence"""
            gge_adapters = {
                5: 'GCAATGAAGACTG',
                3: 'GTGTCTTCTAACG'}
            gge_overhangs = {
                1: {5: 'CCTC', 3: 'CATA'},
                2: {5: 'CATA', 3: 'AAAA'},
                3: {5: 'AAAA', 3: 'AGGA'},
                4: {5: 'AGGA', 3: 'AGTG'},
                5: {5: 'AGTG', 3: 'CAGC'},
                6: {5: 'CAGC', 3: 'TGAA'},
                7: {5: 'TGAA', 3: 'ATTA'},
                8: {5: 'ATTA', 3: 'AATC'},
                9: {5: 'AATC', 3: 'CCAG'}}
            # Determine Left Overhang
            if sequence[:4] == gge_overhangs[modules.left][5]:
                left_overhang = ''
            elif sequence[:3] == gge_overhangs[modules.left][5][1:]:
                left_overhang = gge_overhangs[modules.left][5][:1]
            elif sequence[:2] == gge_overhangs[modules.left][5][2:]:
                left_overhang = gge_overhangs[modules.left][5][:2]
            elif sequence[0] == gge_overhangs[modules.left][5][3:]:
                left_overhang = gge_overhangs[modules.left][5][:3]
            else:
                left_overhang = gge_overhangs[modules.left][5]
            # Determine Right Overhang
            if sequence[-4:] == gge_overhangs[modules.right][3]:
                right_overhang = ''
            elif sequence[-3:] == gge_overhangs[modules.right][3][:3]:
                right_overhang = gge_overhangs[modules.right][3][-1]
            elif sequence[-2:] == gge_overhangs[modules.right][3][:2]:
                right_overhang = gge_overhangs[modules.right][3][-2:]
            elif sequence[-1] == gge_overhangs[modules.right][3][:1]:
                right_overhang = gge_overhangs[modules.right][3][-3:]
            else:
                right_overhang = gge_overhangs[modules.right][3]
            # Append Overhangs
            return gge_adapters[5]+left_overhang+sequence+right_overhang+gge_adapters[3]

        return main()