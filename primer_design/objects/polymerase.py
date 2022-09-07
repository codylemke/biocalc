"""
MODULE: Contains the Polymerase class

Classes:
    Polymerase
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'


# CLASSES ---------------------------------------------------------------------
class Polymerase:
    """Representation of a nucleotide"""
    
    def __init__(
        self,
        buffer,
        buffer_concentration,
        dntp_concentration,
        mg_concentration,
        primer_concentration,
        template_concentration,
        polymerase_concentration,
        extension_rate,
        initial_denaturation,
        denature_step,
        anneal_step,
        extension_step,
        cycles,
        final_extension,
        indefinite_hold,
    ):
        """Initialize pcr reaction instance"""
        self.buffer = buffer
        self.buffer_concentration = buffer_concentration
        self.dntp_concentration = dntp_concentration
        self.mg_concentration = mg_concentration
        self.primer_concentration = primer_concentration
        self.template_concentration = template_concentration
        self.polymerase_concentration = polymerase_concentration
        self.extension_rate = extension_rate
        self.initial_denaturation = initial_denaturation
        self.denature_step = denature_step
        self.anneal_step = anneal_step
        self.extension_step = extension_step
        self.cycles = cycles
        self.final_extension = final_extension
        self.indefinite_hold = indefinite_hold
        return None

    def estimate_run_time(self):
        """
        Return a dictionary that specifies the volumes (value) of each
        component (key) in the reaction mix
        """
        reaction_mix = dict()
        return reaction_mix

    def calculate_thermocycler_conditions(self):
        """
        Return a dictionary that describes the details (values) of each step
        (key) in the the thermocycler run
        """
        thermocycler_run = dict()
        return thermocycler_run


    