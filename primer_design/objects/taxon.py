"""
DOCUMENTATION AND MODULES -----------------------------------------------------

"""
# DJANGO
from django.utils.translation import gettext_lazy as _
from django.urls import reverse_lazy
# PANDAS
import pandas as pd
# FIELDS
from django.db.models import (
    ManyToManyField, OneToOneField, CharField, DateTimeField, IntegerField,
    ForeignKey, URLField, PROTECT, SET_NULL
)
from django.contrib.postgres.fields import ArrayField
# MODELS
from .codon_usage import CodonUsage
from .genetic_code import GeneticCode
# MODULES
from bioinformatics.constants import RESIDUES


"""
MODEL -------------------------------------------------------------------------
"""
class Taxon(Entry):
    """Representation of a phylogenetic classification"""
    # META --------------------------------------------------------------------
    class Meta:
        db_table = 'taxon'
        verbose_name = _('Taxon')
        verbose_name_plural = _('Taxons')
        ordering = ('taxid',)
        get_latest_by = ('-updated_on', '-created_on')
    # FIELDS ------------------------------------------------------------------
    entry_ptr = OneToOneField(
        Entry,
        primary_key=False,
        parent_link=True,
        on_delete=PROTECT,
        related_name='entry_taxon',
        db_column='entry_id',
        verbose_name=_('Biolink ID')
    )
    taxid = IntegerField(
        primary_key=True,
        null=False,
        blank=False,
        unique=True,
        editable=False,
        db_index=False,
        db_column='taxid',
        verbose_name=_('Taxon ID'),
        help_text=_('Unique taxon identifier mapped to NCBI taxid')
    )
    ncbi_api_url = URLField(
        max_length=255,
        null=False,
        blank=False,
        unique=True,
        editable=True,
        db_index=False,
        db_column='ncbi_api_url',
        verbose_name=_('NCBI API URL'),
        help_text=_('URL used to retrieve the entry from NCBI')
    )
    ncbi_created_date = DateTimeField(
        auto_now=False,
        auto_now_add=False,
        null=False,
        blank=False,
        unique=False,
        editable=False,
        db_index=False,
        db_column='ncbi_created_date',
        verbose_name=_('NCBI Created Date'),
        help_text=('Date the entry was initially created on NCBI')
    )
    ncbi_updated_date = DateTimeField(
        auto_now=False,
        auto_now_add=False,
        null=False,
        blank=False,
        unique=False,
        editable=True,
        db_index=False,
        db_column='ncbi_updated_date',
        verbose_name=_('NCBI Updated Date'),
        help_text=('Date the entry was last updated on NCBI')
    )
    ncbi_pubdate = DateTimeField(
        auto_now=False,
        auto_now_add=False,
        null=False,
        blank=False,
        unique=False,
        editable=False,
        db_index=False,
        db_column='ncbi_pubdate',
        verbose_name=_('NCBI Publish Date'),
        help_text=('Date the entry was published on NCBI')
    )
    scientific_name = CharField(
        max_length=255,
        null=False,
        blank=False,
        unique=False,
        editable=True,
        db_index=False,
        db_column='scientific_name',
        verbose_name=_('Scientific Name'),
        help_text=_('Scientific name according to NCBI')
    )
    common_name = CharField(
        max_length=255,
        null=True,
        blank=True,
        unique=False,
        editable=True,
        db_index=False,
        db_column='common_name',
        verbose_name=_('Common Name'),
        help_text=_('Common name according to NCBI')
    )
    other_names = ArrayField(
        base_field=CharField(
            max_length=255,
            null=True,
            blank=True,
            unique=False,
            editable=True,
            db_index=False,
        ),
        null=True,
        blank=True,
        default=list,
        unique=False,
        editable=True,
        db_index=False,
        db_column='other_names',
        verbose_name=_('Other Names'),
        help_text=_('Other names refer to this taxon according to NCBI')
    )
    rank = CharField(
        max_length=255,
        null=False,
        blank=False,
        unique=False,
        editable=True,
        db_index=False,
        db_column='rank',
        verbose_name=_('Rank'),
        help_text=_('Phylogenetic rank of the taxon')
    )
    division = CharField(
        max_length=255,
        null=True,
        blank=True,
        unique=False,
        editable=True,
        db_index=False,
        db_column='division',
        verbose_name=_('Division'),
        help_text=_('Phylogenetic division of the taxon')
    )
    genetic_code = ForeignKey(
        GeneticCode,
        on_delete=PROTECT,
        related_name='genetic_code_taxons',
        related_query_name='genetic_code_taxon',
        null=False,
        blank=False,
        unique=False,
        editable=False,
        db_index=False,
        db_column='genetic_code',
        verbose_name=_('Genetic Code'),
        help_text=_('Genetic code according to NCBI')
    )
    codon_usage = ForeignKey(
        CodonUsage,
        on_delete=PROTECT,
        related_name='codon_usage_taxons',
        related_query_name='codon_usage_taxon',
        null=True,
        blank=True,
        unique=False,
        editable=True,
        db_index=False,
        db_column='codon_usage_table',
        verbose_name=_('Codon Usage Table'),
        help_text=_('Codon usage table of the organism according to Kazusa')
    )
    mitochondrial_genetic_code = ForeignKey(
        GeneticCode,
        on_delete=PROTECT,
        related_name='genetic_code_mitochondrias',
        related_query_name='genetic_code_mitochondria',
        null=True,
        blank=False,
        default=None,
        unique=False,
        editable=False,
        db_index=False,
        db_column='mitochondrial_genetic_code',
        verbose_name=_('Mitochondrial Genetic Code'),
        help_text=_('Mitochondrial genetic code according to NCBI')
    )
    plastidial_genetic_code = ForeignKey(
        GeneticCode,
        on_delete=PROTECT,
        related_name='genetic_code_plastids',
        related_query_name='genetic_code_plastid',
        null=True,
        blank=False,
        default=None,
        unique=False,
        editable=False,
        db_index=False,
        db_column='plastidial_genetic_code',
        verbose_name=_('Plastidial Genetic Code'),
        help_text=_('Plastidial genetic code according to NCBI')
    )
    clades = ManyToManyField(
        'self',
        db_table='organism_clades',
        symmetrical=False,
        blank=True,
        unique=False,
        editable=True,
        db_index=False,
        db_column='clades',
        verbose_name=_('Clades'), 
        help_text=_('Clades that this taxon belongs to')
    )
    category = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='category_taxons',
        related_query_name='category_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='category',
        verbose_name=_('Category'), 
        help_text=_('Category that this taxon belongs to')
    )
    superkingdom = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='superkingdom_taxons',
        related_query_name='superkingdom_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='superkingdom',
        verbose_name=_('Superkingdom'), 
        help_text=_('Superkingdom that this taxon belongs to')
    )
    kingdom = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='kingdom_taxon',
        related_query_name='kingdom_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='kingdom',
        verbose_name=_('Kingdom'),
        help_text=_('Kingdom that this taxon belongs to')
    )
    subkingdom = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subkingdom_taxons',
        related_query_name='subkingdom_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subkingdom',
        verbose_name=_('Subkindom'),
        help_text=_('Subkingdom that this taxon belongs to')
    )
    superphylum = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='superphylum_taxons',
        related_query_name='superphylum_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='superphylum',
        verbose_name=_('Superphylum'), 
        help_text=_('Superphylums that this taxon belongs to')
    )
    phylum = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='phylum_taxons',
        related_query_name='phylum_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='phylum',
        verbose_name=_('Phylum'),
        help_text=_('Phylums that this taxon belongs to')
    )
    subphylum = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subphylum_taxons',
        related_query_name='subphylum_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subphylum',
        verbose_name=_('Subphylum'),
        help_text=_('Subphylums that this taxon belongs to')
    )
    superclass = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='superclass_taxons',
        related_query_name='superclass_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='superclass',
        verbose_name=_('Superclass'),
        help_text=_('Superclass that this taxon belongs to')
    )
    taxon_class = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='taxon_class_taxons',
        related_query_name='taxon_class_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='organism_class',
        verbose_name=_('Class'),
        help_text=_('Class that this taxon belongs to')
    )
    subclass = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subclass_taxons',
        related_query_name='subclass_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subclass',
        verbose_name=_('Subclass'), 
        help_text=_('Subclass that this taxon belongs to')
    )
    infraclass = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='infraclass_taxons',
        related_query_name='infraclass_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='infraclass',
        verbose_name=_('Infraclass'), 
        help_text=_('Infraclass that this taxon belongs to')
    )
    cohort = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='cohort_taxons',
        related_query_name='cohort_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='cohort',
        verbose_name=_('Cohort'),
        help_text=_('Cohort that this taxon belongs to')
    )
    superorder = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='superorder_taxons',
        related_query_name='superorder_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='superorder',
        verbose_name=_('Superorder'),
        help_text=_('Superorder that this taxon belongs to')
    )
    order = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='order_taxons',
        related_query_name='order_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='order',
        verbose_name=_('Order'),
        help_text=_('Order that this taxon belongs to')
    )
    suborder = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='suborder_taxons',
        related_query_name='suborder_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='suborder',
        verbose_name=_('Suborder'), 
        help_text=_('Suborder that this taxon belongs to')
    )
    infraorder = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='infraorder_taxons',
        related_query_name='infraorder_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='infraorder',
        verbose_name=_('Infraorder'),
        help_text=_('Infraorder that this taxon belongs to')
    )
    parvorder =ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='parvorder_taxons',
        related_query_name='parvorder_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='parvorder',
        verbose_name=_('Parvorder'),
        help_text=_('Parvorder that this taxon belongs to')
    )
    superfamily = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='superfamily_taxons',
        related_query_name='superfamily_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='superfamily',
        verbose_name=_('Superfamily'),
        help_text=_('Superfamily that this taxon belongs to')
    )
    family = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='family_taxons',
        related_query_name='family_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='family',
        verbose_name=_('Family'),
        help_text=_('Family that this taxon belongs to')
    )
    subfamily = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subfamily_taxons',
        related_query_name='subfamily_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subfamily',
        verbose_name=_('Subfamily'), 
        help_text=_('Subfamily that this taxon belongs to')
    )
    supertribe = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='supertribe_taxons',
        related_query_name='supertribe_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='supertribe',
        verbose_name=_('Supertribe'),
        help_text=_('Supertribe that this taxon belongs to')
    )
    tribe = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='tribe_taxons',
        related_query_name='tribe_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='tribe',
        verbose_name=_('Tribe'),
        help_text=_('Tribe that this taxon belongs to')
    )
    subtribe = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subtribe_taxons',
        related_query_name='subtribe_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subtribe',
        verbose_name=_('Subtribe'),
        help_text=_('Subtribe that this taxon belongs to')
    )
    genus = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='genus_taxons',
        related_query_name='genus_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='genus',
        verbose_name=_('Genus'),
        help_text=_('Genus that this taxon belongs to')
    )
    subgenus = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subgenus_taxons',
        related_query_name='subgenus_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subgenus',
        verbose_name=_('Subgenus'),
        help_text=_('Subgenus that this taxon belongs to')
    )
    species_group = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='species_group_taxons',
        related_query_name='species_group_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='species_group',
        verbose_name=_('Species Group'),
        help_text=_('Species group that this taxon belongs to')
    )
    species_subgroup = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='species_subgroup_taxons',
        related_query_name='species_subgroup_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='species_subgroup',
        verbose_name=_('Species Subgroup'),
        help_text=_('Species subgroup that this taxon belongs to')
    )
    species = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='species_taxons',
        related_query_name='species_taxon',
        null=True,
        blank=True,
        default=None,
        unique=False,
        editable=True,
        db_index=False,
        db_column='species',
        verbose_name=_('Species'),
        help_text=_('Species that this taxon belongs to')
    )
    subspecies = ForeignKey(
        'self',
        on_delete=SET_NULL,
        related_name='subspecies_taxons',
        related_query_name='subspecies_taxon',
        null=True,
        blank=True,
        unique=False,
        editable=True,
        db_index=False,
        db_column='subspecies',
        verbose_name=_('Subspecies'),
        help_text=_('Subspecies that this taxon belongs to')
    )
    # METHODS -----------------------------------------------------------------
    def __str__(self):
        """Return string representation"""
        return str(_(f'{self.taxid}: {self.scientific_name}'))

    def get_absolute_url(self):
        """Return URL"""
        return reverse_lazy('taxonomy:taxon-detail', kwargs={'pk': self.pk})

    def get_or_fetch(self):
        """Get object from database if it exists. Otherwise fetch from online"""
        try: entry = Taxon.objects.get(taxid=self.taxid)
        except Taxon.DoesNotExist:
            from twisted.internet import reactor, defer
            from scrapy.crawler import CrawlerRunner
            from scrapy.utils.project import get_project_settings
            from biocrawler.spiders import KazusaSpider, NCBITaxonomySpider
            settings = get_project_settings()
            runner = CrawlerRunner(settings)
            @defer.inlineCallbacks
            def crawl():
                # yield runner.crawl(KazusaSpider, taxid=self.taxid)
                yield runner.crawl(NCBITaxonomySpider, taxid=self.taxid)
                reactor.stop()
            crawl()
            reactor.run()
            entry = Taxon.objects.get(taxid=self.taxid)
        # finally:
        return entry

    def is_organism(self):
        if self.rank == 'species' or self.rank == 'subspecies':
            return True
        else:
            return False

    def codons(self):
        """Return the codon usage table as a dataframe"""
        # Joins the GeneticCode and CodonUsage objects and creates a dataframe
        genetic_code_df = pd.DataFrame({
            'codon': self.genetic_code.codons,
            'residue': self.genetic_code.residues,
            'start_codon': self.genetic_code.start_codons,
            'stop_codon': self.genetic_code.stop_codons,
        })
        codon_usage_df = pd.DataFrame({
            'codon': self.codon_usage.codons,
            'frequency': self.codon_usage.frequency,
            'total_number': self.codon_usage.total_number,
        })
        codons = pd.merge(genetic_code_df, codon_usage_df, on='codon', how='outer')
        # Calculates the total number of codons for each residue
        total_residue_codons = dict()
        for residue in RESIDUES:
            value = codons[codons['residue'] == residue]['total_number'].sum()
            total_residue_codons[residue] = value
        # Calculates the relative frequency of each codon and adds it as a row to the dataframe
        relative_frequency = list()
        for index, row in codons.iterrows():
            relative_frequency.append(
                round((row['total_number'] / total_residue_codons[row['residue']]) * 100, 2)
            )
        codons.insert(2, 'relative_frequency', relative_frequency)
        return codons

