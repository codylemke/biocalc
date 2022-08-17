"""Protein

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASS -----------------------------------------------------------------------
class Protein:
    """Representation of a protein"""

    standard_residues = {
        'A': ['Ala', 'Alanine', 89],
        'R': ['Arg', 'Arginine', 174],
        'N': ['Asn', 'Asparagine', 132],
        'D': ['Asp', 'Aspartic acid', 133],
        'C': ['Cys', 'Cysteine', 121],
        'E': ['Glu', 'Glutamic acid', 147],
        'Q': ['Gln', 'Glutamine', 146],
        'G': ['Gly', 'Glycine', 75],
        'H': ['His', 'Histidine', 155],
        'I': ['Ile', 'Isoleucine', 133],
        'L': ['Leu', 'Leucine', 131],
        'K': ['Lys', 'Lysine', 146],
        'M': ['Met', 'Methionine', 149],
        'F': ['Phe', 'Phenylalanine', 165],
        'P': ['Pro', 'Proline', 115],
        'S': ['Ser', 'Serine', 105],
        'T': ['Thr', 'Threonine', 119],
        'W': ['Trp', 'Tryptophan', 204],
        'Y': ['Tyr', 'Tyrosine', 181],
        'V': ['Val', 'Valine', 117]
    }

    ambiguous_residues = {
        'B': ['D', 'N'],
        'Z': ['E', 'Q']
    }
    
    def __init__(self, sequence, name='generic_protein'):
        self.sequence = sequence
        self.name = name
        return

    @property
    def sequence(self):
        """Custom validation for sequence attribute"""
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        from .nucleotide import Nucleotide
        seq = (seq.replace(' ', '').replace('\n', '').replace('-', '')
               .replace('_', '').strip().upper())
        if any(char not in Protein.standard_residues for char in seq):
            if any(char not in Protein.ambiguous_residues for char in seq):
                if all(char in Nucleotide.standard_bases
                       or char in Nucleotide.ambiguous_bases
                       for char in seq):
                    raise ValueError(
                        'The sequence entered appears to be a nuceotide sequence.')
                else:
                    raise ValueError(
                        'The sequence entered does not appear to be a biological sequence.')
        self._sequence = seq
        return

    @property
    def length(self):
        """Return the number of residues in the sequence"""
        residues = len(self.sequence)
        return residues

    @property
    def molecular_weight(self):
        """Returns the kDa of the protein sequence"""
        mws = list()
        for residue in self.sequence:
            if residue in Protein.ambiguous_residues:
                ambiguous_mws = list()
                for possibility in Protein.ambiguous_residues[residue]:
                    ambiguous_mws.append(Protein.standard_residues[possibility][2])
                average = sum(ambiguous_mws) / len(Protein.ambiguous_residues[residue])
                mws.append(average)
            mws.append(Protein.standard_residues[residue][2])
        mw = sum(mws) / self.length
        return mw

    def contains_ambiguity(self):
        """
        Return 'True' if the sequence contains an ambigous
        nucleotide character
        """
        return any(char in Protein.ambiguous_residues for char in self.sequence)
