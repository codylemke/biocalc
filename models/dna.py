"""Contains the DNA object

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
import re
# Third Party Packages
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
        self, sequence, name='generic_dna', part_type='generic_dna'):
        """Constructor Function"""
        super().__init__(sequence=sequence, name=name)
        self.type = part_type
        return

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        seq = seq.replace(' ', '').replace('/n', '').replace('-', '').strip().upper()
        for nucleotide in seq:
            if nucleotide == 'U':
                raise AttributeError('DNA sequence contains U nucleotide')
            if nucleotide not in Nucleotide.nucleotide_codes.keys():
                # if nucleotide in Protein.residue_codes
                    # raise AttributeError('Sequence appears to be a protein sequence')
                raise AttributeError('Sequence contains not nucleotide chatacters')
        self._sequence = seq
        return

    @property
    def complement(self):
        """
        Return the sequence of the complement sequence in the 3' to 5'
        orientation
        """
        complement_sequence = list()
        for index, nucleotide in enumerate(self.sequence, 1):
            if nucleotide == 'A':
                complement_sequence.append('T')
            elif nucleotide == 'T':
                complement_sequence.append('A')
            elif nucleotide == 'G':
                complement_sequence.append('C')
            elif nucleotide == 'C':
                complement_sequence.append('G')
            else:
                raise ValueError(
                    'An invalid nucleotide is present in the sequence at'
                    f'position {index}')
        return ''.join(complement_sequence)

    @property
    def reverse_complement(self):
        """Return the reverse complement sequence"""
        rev_comp = ''.join(self.complement.split().reverse())
        return rev_comp

    def flip_sequence(self):
        """Represents the sequence as the reverse complement"""
        self.sequence = self.reverse_complement
        return

    @staticmethod
    def calculate_nn_tm(
        seq, check=True, strict=False, c_seq=None, shift=0,
        dnac1=25, dnac2=25, selfcomp=False, Na=50, K=0, Tris=0, Mg=1.5,
        dNTPs=0.6, saltcorr=7):
        """
        Calculates melting temperature using BioPython's most up to date
        nearest-neighbor methods.
        """
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
        Calculates melting temperature using BioPython's most up to date
        nearest-neighbor methods.
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
        Calculates melting temperature using BioPython's most up to date
        nearest-neighbor methods.
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

    def pcr_amplify(self, primer_1=0, primer_2=None, project='generic'):
        return # Fragement

    def design_amplification_primers(self, index_1=0, index_2=None, project='generic', left_extension=None, right_extension=None):
        return # (f_primner, r_primer)

    def find_cutters(self, enzymes):
        return # [[index, recognition_site, cut_site]

    def find_golden_gate_overhangs(self, overhangs):
        return # [index, overhang, direction]

    def make_rna(self, strand, index_1=None, index_2=None):
        return # RNA

    def find_orfs(self):
        return # [index, cds]

    def anneal_oligo(self, primer):
        """Finds the template index based on the sequence"""
        # Find if the last 13 bases of the primer anneal to the template
        three_prime_seed_length = 13
        try:
            binding_sites = [index.start() for index in re.finditer(primer.sequence[-three_prime_seed_length:], self.sequence)]
        except ValueError as err:
            raise err('Annealing site not found on template sequence.')
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

    def simulate_agarose_gel(self):
        return # Image


    def align(self, *DNAs):
        return # Alignment
