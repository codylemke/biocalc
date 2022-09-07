# MODULE ----------------------------------------------------------------------
"""
Contains the FragmentOrder object.

Classes:
    PrimerOrder
"""
__author__ = 'Cody Lemke'
__version__ = '2.0.0'

# Standard Library
import re
from pathlib import Path
from sys import ps1
# Third Party Packages
import openpyxl  # https://openpyxl.readthedocs.io/en/stable/
import python_codon_tables as pct  # https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables
# Local Modules
from primer_design.objects import (
    Variant,
    DNA,
    Primer,
    PCRReaction,
    Fragment,
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
class FragmentOrder:
    """Representation of a fragment order"""
    
    def __init__(
        self,
        order_name,
        desired_sequences,
        fragment_templates,
        destination_taxon,
        polymerase,
        purpose,
        method,
    ):
        """Fragment order constructor function"""
        self.order_name = order_name
        self.pcr_reactions = self.generate_pcr_reactions()
        self.restriction_digests = generate_restriction_digests()
        self.desired_sequences = desired_sequences
        self.fragment_templates = fragment_templates
        self.destination_taxon = destination_taxon
        self.polymerase = polymerase
        self.purpose = purpose
        self.method = method

    def parse_desired_sequences(self):
        """Placeholder"""
        

    def generate_fragments(self):
        """Placeholder"""
        for desired_fragment in self.desired_fragments:
            if self.purpose == 'amplification':
                pass
            elif self.purpose == 'assembly':
                if self.method == 'gibson':
                    pass
                elif self.method == 'golden_gate':
                    pass
                elif self.method == 'cpec':
                    pass
                elif self.method == 'hifi':
                    pass
                elif self.method == 'homologous_recombination':
                    pass
                elif self.method == 'infusion':
                    pass
                elif self.method == 'gateway':
                    pass
                elif self.method == 'ta':
                    pass
                elif self.method == 'topo':
                    pass
                elif self.method == 'ligation':
                    pass
                elif self.method == 'linear_ligation':
                    pass
        
        fragment = Fragment(
            name='',
            sequence='',
            purpose='',
            method='',
            source='',
            primers='',
            project='',
        )

    def generate_pcr_reactions(self):
        """Return list of pcr reactions from the order"""
        pcr_reactions = list()
        for variant_code in self.variant_codes:
            aa_position = re.search(r'\d+', variant_code).group()
            wt_residue = variant_code.split(aa_position)[0]
            variation = variant_code.split(aa_position)[1]
            codon_index = (int(aa_position)*3 - 3, int(aa_position)*3)
            # Trick Primers
            if variation == '#':
                reaction_type = 'trick'
                primer_codons = ['NDT','VMA','ATG','TGG']
                possible_codons = list()
                possible_residues = list()
                for codon in primer_codons:
                    possible_codons.extend(parse_degenerate_codon(codon))
                # Create Variants
                for possible_codon in possible_codons:
                    translated_codon = translate_codon(possible_codon)
                    if translated_codon not in possible_residues:
                        possible_residues.append(translated_codon)
                possible_variants = list()
                for possible_residue in possible_residues:
                    possible_variant = Variant(self.enzyme, f'{wt_residue}{aa_position}{possible_residue}')
                    possible_variants.append(possible_variant)
                # Create Primers
                primer_pairs = list()
                for primer_codon in primer_codons:
                    primer_pair = self.design_primers(codon_index=codon_index, codon=primer_codon)
                    primer_pairs.append(primer_pair)
                primers = list()
                for primer_pair in primer_pairs:
                    for primer in primer_pair:
                        primers.append(primer)
            # Smart Primers
            elif variation == '!':
                reaction_type = 'smart'
                primer_codons = ['NDT','VHG','TGG']
                possible_codons = list()
                possible_residues = list()
                for codon in primer_codons:
                    possible_codons.extend(parse_degenerate_codon(codon))
                # Create Variants
                for possible_codon in possible_codons:
                    translated_codon = translate_codon(possible_codon)
                    if translated_codon not in possible_residues:
                        possible_residues.append(translated_codon)
                variants = list()
                for possible_residue in possible_residues:
                    variant = Variant(self.enzyme, f'{wt_residue}{aa_position}{possible_residue}')
                    variants.append(variant)
                # Create Primers
                primer_pairs = list()
                for primer_codon in primer_codons:
                    primer_pair = self.design_primers(codon_index=codon_index, codon=primer_codon)
                    primer_pairs.append(primer_pair)
                primers = list()
                for primer_pair in primer_pairs:
                    for primer in primer_pair:
                        primers.append(primer)
            # Degenerate Primers
            elif len(variation) > 1:
                reaction_type = 'degenerate'
                primer_codon = variation.replace('(', '').replace(')', '')
                possible_codons = parse_degenerate_codon(primer_codon)
                # Create Variants
                possible_residues = list()
                for possible_codon in possible_codons:
                    translated_codon = translate_codon(possible_codon)
                    if translated_codon not in possible_residues:
                        possible_residues.append(translated_codon)
                variants = list()
                for possible_residue in possible_residues:
                    variant = Variant(self.enzyme, f'{wt_residue}{aa_position}{possible_residue}')
                    variants.append(variant)
                # Create Primers
                primer_pair = self.design_primers(codon_index=codon_index, codon=primer_codon)
                primers = list()
                for primer in primer_pair:
                    primers.append(primer)
            # Standard Primers
            else:
                reaction_type = 'standard'
                variant = Variant(self.enzyme, variant_code)
                # Determine best codon
                best_codon_frequency = 0
                for codon in self.codon_table[variation].keys():
                    codon_frequency = self.codon_table[variation][codon]
                    if codon_frequency > best_codon_frequency:
                        best_codon = codon
                        codon_frequency = codon_frequency
                # Create Primers
                primer_pair = self.design_primers(codon_index=codon_index, codon=best_codon)
                primers = list()
                for primer in primer_pair:
                    primers.append(primer)
            # Create PCR Reaction Object
            pcr_reaction = PCRReaction(
                reaction_volume=self.reaction_volume,
                polymerase=self.polymerase,
                template=self.template,
                template_concentration=self.template_concentration,
                primers=primers,
                reaction_type=reaction_type,
            )
            pcr_reactions.append(pcr_reaction)

    def generate_restriction_digests(self):
        pass