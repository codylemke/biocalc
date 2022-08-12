"""Construct

Description.

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
# Local Modules
from .dna import DNA
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class CDS(DNA):
    """Representation of a construct"""

    def __init__(
        self, sequence, *parts, name='generic_construct',
        project='generic_project', cloning_method='unspecified'):
        super.__init__(sequence=sequence, parts=parts, name='generic_construct',
        project='generic_project', cloning_method='unspecified'):
        return