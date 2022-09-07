"""
DOCUMENTATION AND MODULES -----------------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'
# THIRD PARTY PACKAGES
from Bio.SeqUtils import MeltingTemp
# APPLICATION MODULES
from .nucleotide import Nucleotide


"""
CLASS -------------------------------------------------------------------------
"""
class DNA(Nucleotide):
    """Representation of DNA"""
    def __init__(self, sequence, name='generic_dna'):
        super().__init__(sequence, name=name)
        # self.name
        # self.sequence
        # self.length
        # self.gc_content
        return
    
    def complement(self):
        """Return the sequence of the complement sequence in the 3' to 5' orientation"""
        complement_sequence = str()
        for index, nucleotide in enumerate(self.sequence, 1):
            if nucleotide == 'A':
                complement_sequence += 'T'
            elif nucleotide == 'T':
                complement_sequence += 'A'
            elif nucleotide == 'G':
                complement_sequence += 'C'
            elif nucleotide == 'C':
                complement_sequence += 'G'
            else:
                raise ValueError(f'An invalid nucleotide is present in the sequence at position {index}')
        return complement_sequence

    def reverse_complement(self):
        """Return the reverse complement sequence"""
        complement = self.complement()
        reverse_complement = complement[::-1]
        return reverse_complement

    def calculate_nn_tm(
        self,
        template_sequence,
        primer_concentration=25,
        template_concentration=25,
        sodium_concentration=50,
        potassium_concentration=0,
        tris_concentration=0,
        magnesium_concentration=1.5,
        dntp_concentration=0.6,
    ):
        """
        Calculates melting temperature using BioPython's most up to date nearest-neighbor methods.
        """
        tm = int(MeltingTemp.Tm_NN(
            self.sequence,
            check=True,
            strict=False, # Note: thermodynamic data is not available for all mismatches which requires the strict parameter to be set to False to prevent failing in these scenarios
            c_seq=template_sequence,
            shift=0,
            nn_table=MeltingTemp.DNA_NN4, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            tmm_table=MeltingTemp.DNA_TMM1, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            imm_table=MeltingTemp.DNA_IMM1, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            de_table=MeltingTemp.DNA_DE1, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
            dnac1=primer_concentration,
            dnac2=template_concentration,
            selfcomp=False,
            Na=sodium_concentration,
            K=potassium_concentration,
            Tris=tris_concentration,
            Mg=magnesium_concentration,
            dNTPs=dntp_concentration,
            saltcorr=7, # https://biopython.org/docs/latest/api/Bio.SeqUtils.MeltingTemp.html#Bio.SeqUtils.MeltingTemp.Tm_NN
        ))
        return tm

    def save(self):
        """Updates all dependent parameters"""
        super().save()
        return None

