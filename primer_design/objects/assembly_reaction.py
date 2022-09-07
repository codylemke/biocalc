"""
DOCUMENTATION AND MODULES -----------------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

"""
CLASS -------------------------------------------------------------------------
"""
class AssemblyReaction:
    """Representation of an assembly reaction"""
    
    def __init__(
        self,
        type,
        vector_concentration,
        reaction_volume,
        parts,
        polymerase=None,
        restriction_enzyme=None,
        ligase=None,
        kit=None,
    ):
        """Initialization of pcr reaction instances"""
        self.type = type
        self.vector_concentration = vector_concentration
        self.reaction_volume = reaction_volume
        self.parts = parts
        self.polymerase = polymerase
        self.restriction_enzyme = restriction_enzyme
        self.ligase = ligase
        self.kit = kit
        return None

    def generate_construct(self):
        """Return construct produced by an assembly reaction"""
        construct = None
        return construct


    