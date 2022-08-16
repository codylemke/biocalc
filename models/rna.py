"""RNA

Description.

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
class RNA(Nucleotide):
    """Representation of an RNA molecule"""

    def __init__(
        self, sequence, *parts, name='generic_construct',
        project='generic_project', cloning_method='unspecified'):
        super().__init__(sequence=sequence, name=name)
        self.parts = parts
        self.project = project
        self.cloning_method = cloning_method
        return
