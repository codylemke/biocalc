"""Contains the ORF representation

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
# Local Modules
from .nucleotide import Nucleotide
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class ORF(Nucleotide):
    """Representation of an Open Reading Frame"""
    
    def __init__(self, sequence, name='generic_dna'):
        """Constructor Function"""
        super().__init__(sequence=sequence, name=name)
        
        return

    def translate(self):
        """Translates the ORF into its corresponding protein sequence"""
        # Need Translation Tables
        return