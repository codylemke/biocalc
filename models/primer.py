"""Contains the Primer Object

Description.

Usage Example:
    ...example
"""
# Standard Libraries
import pathlib
import re
# Third Party Packages
from Bio.SeqUtils import MeltingTemp
# Local Modules
from .dna import DNA
# Global Constants
GG_EXPRESS_DIR = pathlib.Path(__file__).resolve().parent.parent


# CLASS -----------------------------------------------------------------------
class Primer(DNA):
    """Representation of a primer"""
    def __init__(
        self, sequence, name='generic_primer', template='none',
        direction='universal',  project='Generic', purpose='Unspecified',
        method='Unspecified', template_index=None):
        """Primer constructor function"""
        super().__init__(sequence=sequence, name=name)
        self.purpose = purpose
        self.method = method
        self.project = project    

    def anneal(self, template, required_three_prime_match=10, required_tm=40,
        additional_five_prime_match=15):
        """Returns the binding index, mismatch index, and direction of the
        primer on the specified template.
        """
        if template.length > 1_000_000:
            required_three_prime_match = 15
            required_tm = 40
        # Find all primer binding sites
        try:
            binding_sites = [index.start() for index in re.finditer(self.sequence[-required_three_prime_match:], template.sequence)]
        except ValueError as err:
            raise err('Annealing site not found on template sequence.')
        else:
            for binding_site in binding_sites:
                binding_index = [binding_site, binding_site + three_prime_seed_length]
                # Find how much of the 3' End of the primer anneals to the template
                while True:
                    three_prime_seed_length += 1
                    try:
                        binding_index[0] = template.sequence.index(self.sequence[-three_prime_seed_length:])
                    except ValueError:
                        three_prime_seed_length -= 1
                        annealed_three_prime_sequence = template.sequence[-three_prime_seed_length:]
                        primer_mismatch_index = [0, self.sequence.index(annealed_three_prime_sequence)]
                        break
                    else:
                        if len(template.sequence[-three_prime_seed_length:]) == len(self.sequence):
                            primer_mismatch_index = None
                            return binding_index, primer_mismatch_index
                # Find if the 5' end of the primer binds to the template
                five_prime_seed_length = 13
                try:
                    template.sequence.index(self.sequence[:five_prime_seed_length])
                except ValueError:
                    return binding_index, primer_mismatch_index
                else:
                # Find how much of the 5' end of the primer anneals to the template
                    while True:
                        five_prime_seed_length += 1
                        try:
                            template.sequence.index(self.sequence[:five_prime_seed_length])
                        except ValueError:
                            five_prime_seed_length -= 1
                            annealed_five_prime_sequence = self.sequence[:five_prime_seed_length]
                            binding_index[0] = template.sequence.index(self.sequence[:five_prime_seed_length])
                            primer_mismatch_index[0] = binding_index[0] + len(annealed_five_prime_sequence)
                            return binding_index, primer_mismatch_index


class AnnealedPrimer:
    """Representation of one nucleotide bound to another"""

    def __init__(self):
        return