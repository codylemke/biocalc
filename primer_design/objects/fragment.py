# MODULE ----------------------------------------------------------------------
"""
Contains the Fragment object.

Classes:
    FragmentOrder
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

# Standard Library
import re
from pathlib import Path
# Third Party Packages
import openpyxl  # https://openpyxl.readthedocs.io/en/stable/
import python_codon_tables as pct  # https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables
# Local Modules
from primer_design.objects import (
    Variant,
    DNA,
    Primer,
    PCRReaction,
)
from primer_design.algorithms import (
    design_amplification_primers,
    design_sequencing_primers,
    design_overlapping_primers,
    design_end_to_end_primers,
    design_golden_mutagenesis_primers,
    design_gibson_primers,
    design_cpec_primers,
    design_golden_gate_primers,
    design_hifi_primers
)

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
EXPERIMENTS_DIR = PROJECT_DIR / 'Experiments'


# CLASSES ---------------------------------------------------------------------
class Fragment(DNA):
    """Representation of a fragment order"""
    
    def __init__(self, name, sequence, purpose, method, source, primers, project):
        """Placeholder"""
        super().__init__(sequence, name=name)
        # self.name
        # self.sequence
        # self.length
        # self.gc_content
        self.purpose = purpose
        self.method = method
        self.source = source
        self.primers = primers
        self.project = project

