"""Contains the FASTA object

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
# Local Modules
from .nucleotide import Nucleotide
from .protein import Protein
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Fasta:
    """Representation of a FASTA file."""
    
    def __init__(self, *paths, name='generic_dna', dna_type='generic'):
        """Constructor Function"""
        super().__init__(sequence=sequence, name=name, nucleotide_type='dna')
        self.dna_type = dna_type
        return

    @classmethod
    def combine_fasta_files(cls):
        """Combines the sequences of the fasta files provided into one"""