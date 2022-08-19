"""Contains the PCR Reaction representation

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
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class PCRReaction:
    """Representation of a PCR Reaction"""
    
    def __init__(self, sequence, name='generic_dna'):
        """Constructor Function"""
        super().__init__(sequence=sequence, name=name)
        
        return