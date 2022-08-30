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
    # MWs obtained from pubchem
    # residue, abbreviation, name, daltons
        'A': ['Ala', 'Alanine', 89.09],
        'R': ['Arg', 'Arginine', 174.20],
        'N': ['Asn', 'Asparagine', 132.12],
        'D': ['Asp', 'Aspartic acid', 133.10],
        'C': ['Cys', 'Cysteine', 121.16],
        'E': ['Glu', 'Glutamic acid', 147.13],
        'Q': ['Gln', 'Glutamine', 146.14],
        'G': ['Gly', 'Glycine', 75.07],
        'H': ['His', 'Histidine', 155.15],
        'I': ['Ile', 'Isoleucine', 131.17],
        'L': ['Leu', 'Leucine', 131.17],
        'K': ['Lys', 'Lysine', 146.19],
        'M': ['Met', 'Methionine', 149.21],
        'F': ['Phe', 'Phenylalanine', 165.19],
        'P': ['Pro', 'Proline', 115.13],
        'S': ['Ser', 'Serine', 105.09],
        'T': ['Thr', 'Threonine', 119.12],
        'W': ['Trp', 'Tryptophan', 204.22],
        'Y': ['Tyr', 'Tyrosine', 181.19],
        'V': ['Val', 'Valine', 117.15],
        '*': ['Stop', 'Stop', 0]}
    nonstandard_residues = {
        'U': ['Sec', 'Selenocysteine', 167.06],
        'O': ['Pyl', 'Pyrrolysine', 255.31]}
    ambiguous_residues = {
        'B': ['Asx', ['D', 'N'], 132.61],
        'Z': ['Glx', ['E', 'Q'], 146.64],
        'J': ['Xle', ['L', 'I'], 131.17],
        'X': ['Xaa', ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F',
                      'P','S','T','W','Y','V'], 136.90]}
    
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
    def molecular_weight(self):
        """Returns the molecular weight in kDa of the protein"""
        mws = list()
        for residue in self.sequence:
            if residue in self.ambiguous_residues:
                ambiguous_mws = list()
                for possibility in self.ambiguous_residues[residue]:
                    ambiguous_mws.append(self.standard_residues[possibility][2])
                average = sum(ambiguous_mws) / len(self.ambiguous_residues[residue])
                mws.append(average)
            mws.append(self.standard_residues[residue][2])
        mw = sum(mws) / self.length
        return mw

    def contains_ambiguity(self):
        """Return `True` if the sequence contains an ambigous residue"""
        return any(residue in self.sequence for residue in self.ambiguous_residues)
    
    def contains_nonstandard_residue(self):
        """Return `True` if the sequence contains a nonstandard residue"""
        return any(residue in self.sequence for residue in self.nonstandard_residues)
