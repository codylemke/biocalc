"""Contains the Nucleotide Object

Description...

Usage Example:
    ...example
"""
# Standard Library
import pathlib
import textwrap
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Nucleotide:
    """Representation of a nucleotide"""

    standard_bases = {
        'A': ['Adenine', 'Adenosine', 0],
        'C': ['Cytosine', 'Cytidine', 0],
        'G': ['Guanine', 'Guanosine', 0],
        'T': ['Thymine', 'Thymidine', 0],
        'U': ['Uracil', 'Uridine', 0]}
    ambiguous_bases = {
    'W': ['A', 'T'],
    'S': ['C', 'G'],
    'M': ['A', 'C'],
    'K': ['G', 'T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']}

    # def __new__(cls, *args, **kwargs):
    #     """Figure out how to auto generate names"""
    #     instance = super().__new__(cls)
    #     return instance


    def __init__(self, sequence, name='generic_nucleotide', nucleotide_type='generic'):
        self.name = name
        self.nucleotide_type = nucleotide_type
        self.sequence = sequence
        return

    @property
    def sequence(self):
        """Requires docstring"""
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        from .protein import Protein
        seq = (seq.replace(' ', '').replace('\n', '').replace('-', '')
               .replace('_', '').strip().upper())
        if not all(char in Nucleotide.standard_bases
                   or char in Nucleotide.ambiguous_bases
                   for char in seq):
            if all(char in Protein.standard_residues
                   or char in Protein.ambiguous_residues
                   for char in seq):
                raise ValueError(
                    'The sequence entered appears to be a protein sequence.')
            else:
                raise ValueError(
                    'The sequence entered does not appear to be biological.')
        if self.nucleotide_type == 'generic':
            if 'U' in seq and 'T' not in seq:
                self.nucleotide_type = 'RNA'
            elif 'T' in seq and 'U' not in seq:
                self.nucleotide_type = 'DNA'
            else:
                raise AttributeError(textwrap.dedent("""\
                    The `nucleotide_type` could not be determined.
                    Please manually specify the `nucleotide_type` as a keyword argument.
                    eg. DNA or RNA"""))
        elif self.nucleotide_type == 'DNA':
            if 'U' in seq and 'T' not in seq:
                raise AttributeError(
                    'The `nucleotide_type` attribute was '
                    'specified as DNA but appears to be RNA')
        elif self.nucleotide_type == 'RNA':
            if 'T' in seq and 'U' not in seq:
                raise AttributeError(
                    'The `nucleotide_type` attribute was '
                    'specified as RNA but appears to be DNA')
        self._sequence = seq
        return

    @property
    def length(self):
        """Return the number of base pairs in the sequence"""
        bases = len(self.sequence)
        return bases

    @property
    def molecular_weight(self):
        """Returns the molecular weight based on sequence"""
        mw = 0
        return mw

    @property
    def gc_content(self):
        """Return the GC content of the sequence"""
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        at_count = self.sequence.count('A') + self.sequence.count('T') + self.sequence.count('U')
        if any(char in Nucleotide.ambiguous_bases for char in self.sequence):
            for ambiguity, chance in Nucleotide.ambiguous_bases.items():
                ambiguity_count = self.sequence.count(ambiguity)
                gc_chance = (chance.count('G') + chance.count('C')) / len(chance)
                gc_count += gc_chance * ambiguity_count
                at_count += ambiguity_count - gc_chance
        percent = round(float(gc_count / (gc_count + at_count) * 100), 2)
        return percent


    def reverse(self):
        """Return the sequence in the 3' to 5' orientation"""
        reverse_sequence = ''.join(reversed(list(self.sequence)))
        return reverse_sequence

    @property
    def complement(self):
        """
        Return the sequence of the complement sequence in the 3' to 5'
        orientation
        """
        complement_sequence = list()
        for nucleotide in self.sequence:
            if nucleotide == 'A':
                if self.nucleotide_type == 'DNA':
                    complement_sequence.append('T')
                elif self.nucleotide_type == 'RNA':
                    complement_sequence.append('U')
            elif nucleotide == 'T' or nucleotide == 'U':
                complement_sequence.append('A')
            elif nucleotide == 'G':
                complement_sequence.append('C')
            elif nucleotide == 'C':
                complement_sequence.append('G')
            elif nucleotide == 'W':
                complement_sequence.append('W')
            elif nucleotide == 'S':
                complement_sequence.append('S')
            elif nucleotide == 'M':
                complement_sequence.append('K')
            elif nucleotide == 'K':
                complement_sequence.append('M')
            elif nucleotide == 'R':
                complement_sequence.append('Y')
            elif nucleotide == 'Y':
                complement_sequence.append('R')
            elif nucleotide == 'B':
                complement_sequence.append('V')
            elif nucleotide == 'D':
                complement_sequence.append('H')
            elif nucleotide == 'H':
                complement_sequence.append('D')
            elif nucleotide == 'V':
                complement_sequence.append('B')
            elif nucleotide == 'N':
                complement_sequence.append('N')
            else:
                raise ValueError(
                    'An invalid nucleotide is present in the sequence')
        return ''.join(complement_sequence)

    @property
    def reverse_complement(self):
        """Return the reverse complement sequence"""
        rev_comp = ''.join(reversed(list(self.complement)))
        return rev_comp

    def flip_sequence(self):
        """Represents the sequence as the reverse complement"""
        self.sequence = self.reverse_complement
        return

    def contains_ambiguity(self):
        """
        Return 'True' if the sequence contains an ambigous
        nucleotide character
        """
        return any(char in Nucleotide.ambiguous_bases for char in self.sequence)
