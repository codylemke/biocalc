"""Contains the Nucleotide Object

Description...

Usage Example:
    ...example
"""
# Standard Library
import pathlib
# Local Modules
from .sequence import Sequence
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Nucleotide(Sequence):
    """Representation of a nucleotide"""

    standard_bases = {
    #  Base: [Nucleotide, Nucleoside, RNMP mw]
        'A': ['Adenine', 'Adenosine', 347.22],
        'C': ['Cytosine', 'Cytidine', 323.20],
        'G': ['Guanine', 'Guanosine', 363.22],
        'T': ['Thymine', 'Thymidine', 338.19],
        'U': ['Uracil', 'Uridine', 324.18]}
    nonstandard_bases = {
        '(i)': ['Hypoxanthene', 'Inosine', 348.21]}
    ambiguous_bases = {
        'R': ['A', 'G'], # purine
        'Y': ['C', 'T', 'U'], # pyrimidine
        'W': ['A', 'T', 'U'], # weak interaction
        'S': ['C', 'G'], # strong interaction
        'M': ['A', 'C'], # amino groups
        'K': ['G', 'T', 'U'], # ketones
        'B': ['C', 'G', 'T', 'U'], # not A
        'D': ['A', 'G', 'T', 'U'], # not C
        'H': ['A', 'C', 'T', 'U'], # not G
        'V': ['A', 'C', 'G'], # not T/U
        'N': ['A', 'C', 'G', 'T', 'U']} # any nucleotide

    # def __new__(cls, *args, **kwargs):
    #     """Figure out how to create the right object regardless of what the user selects"""
    #     instance = super().__new__(cls)
    #     return instance
    # recursion can be used when saving inherited models to the DB

    def __init__(self, sequence, name='generic_nucleotide', strands=None,
                 nucleotide_type=None, orientation='linear'):
        self.nucleotide_type = nucleotide_type
        super().__init__(sequence=sequence, name=name, sequence_type='dna')
        self.name = name
        self.sequence = sequence
        self.orientation = orientation
        if strands is None:
            if self.nucleotide_type == 'DNA':
                self.strands = 'double'
            elif self.nucleotide_type == 'RNA':
                self.strands = 'single'
        else:
            self.strands = strands
        return

    @property
    def sequence(self):
        """Placeholder"""
        return Sequence.sequence.fget(self)

    @sequence.setter
    def sequence(self, seq):
        """Placeholder"""
        Sequence.sequence.fset(self, seq)
        if self.nucleotide_type is None:
            if 'U' in self.sequence and 'T' not in self.sequence:
                self.nucleotide_type = 'RNA'
            elif 'T' in self.sequence and 'U' not in self.sequence:
                self.nucleotide_type = 'DNA'
            else:
                raise AttributeError(
                    'The `nucleotide_type` could not be determined')
        elif self.nucleotide_type == 'DNA':
            if 'U' in self.sequence and 'T' not in self.sequence:
                raise AttributeError(
                    'The `nucleotide_type` attribute was '
                    'specified as DNA but appears to be RNA')
        elif self.nucleotide_type == 'RNA':
            if 'T' in self.sequence and 'U' not in self.sequence:
                raise AttributeError(
                    'The `nucleotide_type` attribute was '
                    'specified as RNA but appears to be DNA')
        return

    @property
    def molecular_weight(self):
        """Returns the molecular weight in g/mol of the sequence"""
        mw = int()
        for standard_base, base_details in self.standard_bases.items():
            if self.nucleotide_type == 'DNA':
                mw += (base_details[2] - 34) * self.sequence.count(standard_base)
                if self.orientation == 'linear':
                    mw += 79
                elif self.orientation == 'circular':
                    pass
                else:
                    pass
            elif self.nucleotide_type == 'RNA':
                mw += (base_details[2] - 18) * self.sequence.count(standard_base)
                if self.orientation == 'linear':
                    mw += 159
                elif self.orientation == 'circular':
                    pass
                else:
                    pass
            else:
                pass
        if not all(char in self.standard_bases for char in self.sequence):
            if self.contains_nonstandard_base():
                for nonstandard_base, base_details in self.nonstandard_bases.items():
                    if self.nucleotide_type == 'DNA':
                        mw += (base_details[2] - 34) * self.sequence.count(nonstandard_base)
                    elif self.nucleotide_type == 'RNA':
                        mw += (base_details[2] - 18) * self.sequence.count(nonstandard_base)
            if self.contains_ambiguity():
                for ambiguous_base, possible_bases in self.ambiguous_bases.items():
                    average_mw = sum(self.standard_bases[possible_base][2] for possible_base in possible_bases) / self.length
                    mw += average_mw * self.sequence.count(ambiguous_base)
        if self.strands == 'double':
            mw *= 2
        return mw

    @property
    def gc_content(self):
        """Return the GC content of the sequence"""
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        at_count = self.sequence.count('A') + self.sequence.count('T') + self.sequence.count('U')
        if self.contains_ambiguity():
            for ambiguous_base, possibilities in Nucleotide.ambiguous_bases.items():
                ambiguous_base_count = self.sequence.count(ambiguous_base)
                gc_chance = (possibilities.count('G') + possibilities.count('C')) / len(possibilities)
                gc_count += gc_chance * ambiguous_base_count
                at_count += (1 - gc_chance) * ambiguous_base_count
        percent = round(float(gc_count / (gc_count + at_count) * 100), 2)
        return percent

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
            elif nucleotide == '(i)':
                complement_sequence.append('N')
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
    def reverse(self):
        """Return the sequence in the 3' to 5' orientation"""
        reverse_sequence = ''.join(reversed(list(self.sequence)))
        return reverse_sequence

    @property
    def reverse_complement(self):
        """Return the reverse complement sequence"""
        rev_comp = ''.join(reversed(list(self.complement)))
        return rev_comp

    def contains_ambiguity(self):
        """Return `True` if the sequence contains an ambigous base"""
        return any(base in self.sequence for base in self.ambiguous_bases)

    def contains_nonstandard_base(self):
        """Returns `True` if the sequence contains a nonstandard base"""
        return any(base in self.sequence for base in self.nonstandard_bases)

    def flip_sequence(self):
        """Represents the sequence as the reverse complement"""
        self.sequence = self.reverse_complement
        return

    @staticmethod
    def mafft_align(sequence_objects):
        """Generates a multiple sequence alignment using mafft and returns
        the location of the output file"""
        import subprocess
        fasta_string = Nucleotide.create_fasta(sequence_objects)
        sequences_path = ROOT_DIR / 'static' / 'temp' / 'sequences.fasta'
        with open(sequences_path, mode='w', encoding='utf-8') as sequences_fasta:
            sequences_fasta.write(fasta_string)
        alignment_path = ROOT_DIR / 'static' / 'temp' / 'alignment.fasta'
        alignment_path.touch()
        subprocess.run(f'mafft --auto --maxiterate 100 --thread 4 {str(sequences_path)} > {str(alignment_path)}', check=True, shell=True)
        with open(alignment_path, mode='r', encoding='utf-8') as alignment_fasta:
            alignment_string = alignment_fasta.read()
        # sequences_path.unlink()
        # alignment_path.unlink()
        return alignment_string
