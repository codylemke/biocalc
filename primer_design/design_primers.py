"""
SCRIPT:
-------------------------------------------------------------------------------
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

# STANDARD LIBRARIES
from pathlib import Path
# LOCAL MODULES
from primer_design.algorithms import (
    design_amplification_primers,
    design_sequencing_primers,
    design_overlapping_primers,
    design_end_to_end_primers,
    design_golden_mutagenesis_primers,
    design_gibson_primers,
    design_cpec_primers,
    design_golden_gate_primers,
)
from primer_design.objects import PrimerOrder

# CONFIGURATION
PROJECT_DIR = Path(__file__).resolve().parent.parent
EXPERIMENTS_DIR = PROJECT_DIR / 'Experiments'


# SCRIPT ----------------------------------------------------------------------
def main():
    """Execution of application logic"""
    experiment_name = input('Enter experiment name:  ')
    order = PrimerOrder(experiment_name)
    if order.purpose == 'amplification':
        if order.method == 'amplification':
            design_amplification_primers(order)
        else:
            raise ValueError('Invalid amplification design method')
    elif order.purpose == 'sequencing':
        if order.method == 'sequencing':
            design_sequencing_primers(order)
        else:
            raise ValueError('Invalid sequencing design method')
    elif order.purpose == 'mutagenesis':
        if order.method == 'overlapping':
            design_overlapping_primers(order)
        elif order.method == 'end_to_end':
            design_end_to_end_primers(order)
        elif order.method == 'golden_mutagenesis':
            design_golden_mutagenesis_primers(order)
        else:
            raise ValueError('Invalid mutagenesis design method')
    elif order.purpose == 'assembly':
        if order.method == 'gibson':
            design_gibson_primers(order)
        elif order.method == 'cpec':
            design_cpec_primers(order)
        elif order.method == 'golden_gate':
            design_golden_gate_primers()
        else:
            raise ValueError('Invalid assembly design method')
    
    else:
        raise ValueError('Invalid primer purpose')


# SCRIPT EXECUTION ------------------------------------------------------------
if __name__ == "__main__":
    main()

