"""Protein

Description...

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASS -----------------------------------------------------------------------
class Protein:
    """Representation of a primer"""
    
    residue_codes = [
        ['A', 'Ala', 'Alanine'],
        ['R', 'Arg', 'Arginine'],
        ['N', 'Asn', 'Asparagine'],
        ['D', 'Asp', 'Aspartic acid'],
        ['C', 'Cys', 'Cysteine'],
        ['E', 'Glu', 'Glutamic acid'],
        ['Q', 'Gln', 'Glutamine'],
        ['G', 'Gly', 'Glycine'],
        ['H', 'His', 'Histidine'],
        ['I', 'Ile', 'Isoleucine'],
        ['L', 'Leu', 'Leucine'],
        ['K', 'Lys', 'Lysine'],
        ['M', 'Met', 'Methionine'],
        ['F', 'Phe', 'Phenylalanine'],
        ['P', 'Pro', 'Proline'],
        ['S', 'Ser', 'Serine'],
        ['T', 'Thr', 'Threonine'],
        ['W', 'Trp', 'Tryptophan'],
        ['Y', 'Tyr', 'Tyrosine'],
        ['V', 'Val', 'Valine']
    ]
    
    def __init__():
        return