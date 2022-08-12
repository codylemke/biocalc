"""Construct

Description.

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
# Third Party Packages
from Bio.SeqUtils import MeltingTemp
# Local Modules
from .dna import DNA
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Construct(DNA):
    """Representation of a construct"""

    def __init__(
        self, sequence, *parts, name='generic_construct',
        project='generic_project', cloning_method='unspecified'):
        """Construct constructor"""
        super().__init__(sequence=sequence, name=name)
        self.project = project
        self.cloning_method = cloning_method
        self.map_created=None
        self.designed_by=None
        self.location=None
        self.sequence_verified=None
        self.cloned_by=None
        self.date_banked=None
        self.banked_by=None
        self.strain=None
        self.notes=None
        self.parts=parts
        return None

    # Use properties to alter how functions interact with the sequence

    @classmethod
    def cpec_assemble(cls, *fragments):
        return # Construct

    @classmethod
    def gibson_assemble(cls, *fragments):
        
        return # Construct

    @classmethod
    def hifi_assemble(cls, *fragments):
        return # Construct

    @classmethod
    def in_fusion_assemble(cls, *fragments):
        return # Construct

    @classmethod
    def golden_gate_assemble(cls, *dnas, enzyme):
        return # Construct

    def linearize(self, index):
        return # DNA

    def set_origin(self, index):
        
        return # None

    