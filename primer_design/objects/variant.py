"""
DOCUMENTATION AND CONFIGURATION -----------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'
# STANDARD LIBRARIES
import re


"""
OBJECT ------------------------------------------------------------------------
"""
class Variant:
    """Representation of an enzyme variation"""
    def __init__(self, enzyme, variant_code):
        """Variant constructor function"""
        self.enzyme = enzyme
        self.variant_code = variant_code
        self.parse_variant_code(variant_code)
        # self.aa_position
        # self.wt_residue
        # self.variant_residue
        return None

    def parse_variant_code(self, variant_code):
        """Returns parsed variant details from variant code"""
        self.aa_position = int(re.search(r'\d+', variant_code).group())
        self.wt_residue = variant_code.split(str(self.aa_position))[0]
        self.variant_residue = variant_code.split(str(self.aa_position))[1].strip()
        return None

