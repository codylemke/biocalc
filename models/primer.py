"""Contains the Primer Object

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


# CLASS -----------------------------------------------------------------------
class Primer(DNA):
    """Representation of a primer"""
    def __init__(
        self, seq, name='generic_primer', template='none',
        direction='universal',  project='Generic', purpose='Unspecified',
        method='Unspecified', template_index=None):
        """Primer constructor function"""
        super().__init__(seq=seq, name=name)
        self.template = template
        self._template_index = template_index
        self.purpose = purpose
        self.method = method
        self.project = project
        self.direction = direction      

    @property
    def template_index(self):
        return self._template_index

    @template_index.setter
    def template_index(self, index):
        self._index = index
        if self.direction == 'F':
            self.sequence = self.template[index[0]:index[1]]
        elif self.direction == 'R':
            self.sequence = self.complement[index[0]:index[1]:-1]

    @property
    def tm(self):
        tm_i = self.calculate_nn_tm(self.sequence, self.template.complement, shift=self.template_index[0])
        return tm_i

    @property
    def tm2(self):
        tm_f = self.calculate_nn_tm(self.sequence, self.sequence.complement + self.template[self.index[1]])
        return tm_f
