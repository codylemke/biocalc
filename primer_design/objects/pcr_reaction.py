"""
MODULE: Contains the PCRReaction class

Classes:
    PCRReaction
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

# Local Modules
from polymerase import Polymerase
from 


# CLASSES ---------------------------------------------------------------------
class PCRReaction:
    """Representation of a nucleotide"""
    
    def __init__(
        self,
        reaction_volume,
        polymerase,
        template,
        template_concentration,
        primers,
        reaction_type,
    ):
        """Initialize pcr reaction instance"""
        self.reaction_volume = reaction_volume
        self.polymerase = polymerase
        self.template = template
        self.template_concentration = template_concentration
        self.template_volumne = self.calculate_template_volume()
        self.primers = primers
        self.reaction_type = reaction_type
        self.reaction_mix = self.calculate_reaction_mix()
        self.thermocycler_conditions = self.calculate_thermocycler_conditions()
        self.stock_dntp_concentration = ''
        self.stock_mg_concentration = ''
        self.gc_enhancer_volume = ''
        self.dmso_volume = ''
        return None

    def calculate_reaction_mix(self):
        """
        Return a dictionary that specifies the volumes (value) of each
        component (key) in the reaction mix
        """
        buffer_volume = self.reaction_volume / self.polymerase_buffer_concentration
        dntp_volume = self.polymerase.dntp_concentration * self.reaction_volume / self.stock_dntp_concentration
        mg_volume = self.polymerase.mg_concentration * self.reaction_volume / self.stock_mg_concentration
        polymerase_volume = self.polymerase.concentration * self.reaction_volume / self.stock_polymerase_concentration
        primers_volume = self.polymerase.primer_concentration * self.reaction_volume / self.stock_primers_volume
        water_volume = self.reaction_volume - (buffer_volume + dntp_volume + mg_volume + polymerase_volume + primers_volume)
        reaction_mix = dict()
        return reaction_mix

    def calculate_extension_duration(self):
        """
        Return a dictionary that describes the details (values) of each step
        (key) in the the thermocycler run
        """
        extension_duration = self.template.length * self.polymerase.rate
        return extension_duration

    def calculate_template_moles(self):
        """
        Return the volume of template that will be added to the reaction based on the
        concentration of the purified template.
        """
        # https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
        a_mw = 331.2
        t_mw = 322.2
        c_mw = 307.2
        g_mw = 347.2
        at_count = self.template.count('A') + self.template.count('T')
        gc_count = self.template.count('G') + self.template.count('C')
        template_mw = (
            a_mw * at_count
            + t_mw * at_count
            + c_mw * gc_count
            + g_mw * gc_count
        )
        template_moles = self.template_concentration / template_mw
        return template_moles

    