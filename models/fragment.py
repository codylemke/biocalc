"""Contains the Fragment object.

Classes:
    FragmentOrder
"""
# Standard Library
import pathlib
# Local Modules
from .dna import DNA
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASSES ---------------------------------------------------------------------
class Fragment(DNA):
    """Representation of a fragment order"""

    def __init__(
        self, sequence, name, purpose=None, method=None, source=None,
        primers=None, project=None, fragment_type=None):
        """Placeholder"""
        super().__init__(sequence=sequence, name=name, dna_type='fragment')
        self.purpose = purpose
        self.method = method
        self.source = source
        self.primers = primers
        self.project = project
        self.fragment_type = fragment_type

    def ligate(self, *fragments):
        return # DNA
