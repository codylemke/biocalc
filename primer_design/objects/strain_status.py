# MODULE ----------------------------------------------------------------------
"""
Contains StrainStatus Object

Classes:
    StrainStatus
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'
# Local Modules
from .dna import DNA


# CLASS -----------------------------------------------------------------------
class StrainStatus:
    """Representation of a strain status"""

    def __init__(self, strain):
        """Initializer"""
        self.strain = strain
        self.id = id
        self.transformation_successful = False
        self.