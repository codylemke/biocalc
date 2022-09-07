"""
DOCUMENTATION AND CONFIGURATION -----------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'
# APPLICATION MODULES
from .dna import DNA


"""
OBJECT ------------------------------------------------------------------------
"""
class Primer(DNA):
    """Representation of a primer"""
    def __init__(self, sequence, project, purpose, method, index, mismatch_index, template, name='generic_primer', order=None):
        """Primer constructor function"""
        super().__init__(sequence, name=name)
        # self.name
        # self.sequence
        # self.length
        # self.gc_content
        self.project = project
        self.purpose = purpose
        self.method = method
        self.index = index
        self.mismatch_index = mismatch_index
        self.template = template
        self.save()
        # self.tm
        # self.tm2
        # self.five_prime_edge_length
        # self.three_prime_edge_length
        # self.five_prime_edge_tm
        # self.three_prime_edge_tm

    def save(self):
        """
        Updates all dependent parameters upon alteration of independent parameters
        """
        super().save()
        self.tm = self.calculate_nn_tm(self.sequence, self.template.complement()[self.index[0]-1:self.index[1]+1], shift=1)
        self.tm2 = self.calculate_nn_tm(self.sequence, f'{self.template[self.index[0]-1]}{self.sequence.complement()}{self.template[self.index[1]]}', shift=1)
        self.five_prime_edge_length = len(self.template[self.index[0]:self.mismatch_index[0]])
        self.three_prime_edge_length = len(self.template[self.mismatch_index[1]:self.index[1]])
        self.five_prime_edge_tm = self.calculate_nn_tm(self.template[self.index[0]:self.mismatch_index[0]], self.template.complement()[self.index[0]-1:self.mismatch_index[0]+1], shift=1)
        self.three_prime_edge_tm = self.calculate_nn_tm(self.template[self.mismatch_index[1]:self.index[1]], self.template.complement()[self.mismatch_index[1]-1:self.index[1]+1], shift=1)
        return None