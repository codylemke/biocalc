"""Contains the Nucleotide Object

Description...

Usage Example:
    ...example
"""
# Standard Library
import pathlib
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Nucleotide:
    """Representation of a nucleotide"""

    standard_nucleotides = ['A','C','G','T','U']

    ambiguous_nucleotides = {
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

    # Create class constructor which auto names nucleotides

    def __init__(self, sequence, name='generic_nucleotide', nucleotide_type=None):
        self.name = name
        self._sequence = sequence
        if nucleotide_type is None:
            if 'U' in self.sequence and 'T' not in self.sequence:
                self.type = 'RNA'
            elif 'T' in self.sequence and 'U' not in self.sequence:
                self.type = 'DNA'
            else:
                raise AttributeError(
                    'The `nucleotide_type` could not be determined.'
                    'Please manually specify the `nucleotide_type` as a keyword argument.'
                    'eg. DNA or RNA')
        return

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        seq = (seq.replace(' ', '').replace('/n', '').replace('-', '')
            .replace('_', '').strip().upper())
        if any(seq.split()) not in Nucleotide.standard_nucleotides:
            if any(['J', 'O', 'X', 'Z']) in seq:
                raise ValueError(
                    'The sequence entered is not a biological sequence')
            elif any(['E', 'F', 'I', 'L', 'P', 'Q']) in seq:
                if not any(['B', 'U']) in seq:
                    raise ValueError(
                        'The sequence entered appears to be a protein '
                        'sequence.')
                else:
                    raise ValueError(
                        'The sequence entered contains non-nucleotide '
                        'characters')
        else:
            self._sequence = seq
        return

    @property
    def length(self):
        """Return the number of base pairs in the sequence"""
        bp = len(self._sequence)
        return bp

    @property
    def gc_content(self):
        """Return the GC content of the sequence"""
        # make this able to handle ambiguity
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        at_count = self.sequence.count('A') + self.sequence.count('T') + self.sequence.count('U')
        if any(Nucleotide.ambiguous_nucleotides) in self.sequence:
            for ambiguity in Nucleotide.ambiguous_nucleotides:
                ambiguity_count = self.sequence.count(ambiguity)
                gc_chance = (
                    (Nucleotide.ambiguous_nucleotides[ambiguity].count('G')
                    + Nucleotide.ambiguous_nucleotides[ambiguity].count('C'))
                    / len(Nucleotide.ambiguous_nucleotides[ambiguity]))
                gc_count += gc_chance * ambiguity_count
                at_count += ambiguity_count - gc_chance
        percent = round(float(gc_count / (gc_count + at_count) * 100), 2)
        return percent


    def reverse(self):
        """Return the sequence in the 3' to 5' orientation"""
        reverse_sequence = ''.join(self.sequence.split().reverse())
        return reverse_sequence


    def contains_ambiguity(self):
        """
        Return 'True' if the sequence contains an ambigous
        nucleotide character
        """
        for ambiguous_nucleotide in Nucleotide.ambiguous_nucleotides:
            if ambiguous_nucleotide in self.sequence:
                return True
        return False
