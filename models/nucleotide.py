"""Contains the Nucleotide Object
"""
# Standard Library
import pathlib
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Nucleotide:
    """Representation of a nucleotide"""
    
    nucleotide_codes = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'U': ['U'],
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
    'N': ['A', 'C', 'G', 'T']
}
    def __init__(self, sequence, name='generic_nucleotide'):
        self.name = name
        self._sequence = sequence
        return

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        seq = seq.strip().upper()
        for nucleotide in seq:
            if nucleotide not in Nucleotide.nucleotide_codes.keys():
                raise ValueError('Sequence contains not nucleotide chatacters')
        self._sequence = seq
        return

    @property
    def length(self):
        """Return the number of base pairs in the sequence"""
        bp = len(self.sequence)
        return bp

    @property
    def gc_content(self):
        """Return the GC content of the sequence"""
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        at_count = self.sequence.count('A') + self.sequence.count('T')
        gc = round(float(gc_count / (gc_count + at_count) * 100), 2)
        return gc

    def reverse(self):
        """Return the sequence in the 3' to 5' orientation"""
        reverse_sequence = ''.join(self.sequence.split().reverse())
        return reverse_sequence
