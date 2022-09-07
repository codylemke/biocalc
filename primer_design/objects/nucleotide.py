"""
DOCUMENTATION AND MODULES -----------------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

"""
CLASS -------------------------------------------------------------------------
"""
class Nucleotide:
    """Representation of a nucleotide"""
    def __init__(self, sequence, name='generic_nucleotide'):
        self.name = name
        self.sequence = sequence
        self.save()
        # self.length
        # self.gc_content
        return

    def reverse(self):
        """Return the sequence in the 3' to 5' orientation"""
        reverse_sequence = self.sequence[::-1]
        return reverse_sequence

    def get_length(self):
        """Return the number of base pairs in the sequence"""
        bp = len(self.sequence)
        return bp

    def get_gc_content(self):
        """Return the GC content of the sequence"""
        gc_count = int()
        at_count = int()
        for index, nucleotide in enumerate(self.sequence, 1):
            if nucleotide == 'G' or nucleotide == 'C':
                gc_count += 1
            elif nucleotide == 'A' or nucleotide == 'T':
                at_count += 1
            else:
                raise ValueError(f'An invalid nucleotide is present in the sequence at position {index}')
        gc_content = round(float(gc_count / (gc_count + at_count) *100), 2)
        return gc_content

    def save(self):
        """Updates all dependent parameters"""
        self.length = self.get_length()
        self.gc_content = self.get_gc_content()
        return