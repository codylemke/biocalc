"""
DOCUMENTATION AND MODULES -----------------------------------------------------

For an overview of this file, see
https://docs.djangoproject.com/en/dev/topics/db/models/

For a comprehensive understanding of this file, see
https://docs.djangoproject.com/en/dev/ref/models/
"""
# DJANGO
from django.utils.translation import gettext_lazy as _
from django.urls import reverse_lazy
# FIELDS
from django.db.models import (
    CharField, IntegerField, BooleanField, OneToOneField, PROTECT,
)
from django.contrib.postgres.fields import ArrayField
# MODELS
from backend.biolink.models import Entry


"""
MODELS ------------------------------------------------------------------------
"""
class GeneticCode(Entry):
    """Representation of a genetic code translation table"""
    # META --------------------------------------------------------------------
    class Meta:
        db_table = 'genetic_code'
        verbose_name = _('Genetic Code')
        verbose_name_plural = _('Genetic Codes')
        ordering = ('genetic_code_id',)
        get_latest_by = ('-updated_on', '-created_on')
    # FIELDS ------------------------------------------------------------------
    entry_ptr = OneToOneField(
        Entry,
        primary_key=False,
        parent_link=True,
        on_delete=PROTECT,
        related_name='entry_genetic_code',
        db_column='entry_id',
    )
    genetic_code_id = IntegerField(
        primary_key=True,
        null=False,
        blank=False,
        unique=True,
        editable=False,
        db_index=False,
        db_column='genetic_code_id',
        verbose_name=_('Genetic Code ID'),
        help_text=_('Unique genetic code identifier')
    )
    name = CharField(
        max_length=255,
        null=False,
        blank=False,
        unique=True,
        editable=True,
        db_index=False,
        db_column='name',
        verbose_name=_('Name'),
        help_text=_('Name of the genetic code according to NCBI')
    )
    codons = ArrayField(
        base_field=CharField(
            max_length=3,
            null=False,
            blank=False,
            unique=True,
            editable=True,
            db_index=False,
        ),
        size=64,
        null=False,
        blank=False,
        default=list,
        unique=False,
        editable=True,
        db_index=False,
        db_column='codons',
        verbose_name=_('Codons'),
        help_text=_('Array of all possible codons')
    )
    residues = ArrayField(
        base_field=CharField(
            max_length=1,
            null=False,
            blank=False,
            unique=False,
            editable=True,
            db_index=False,
        ),
        size=64,
        null=False,
        blank=False,
        default=list,
        unique=False,
        editable=True,
        db_index=False,
        db_column='residues',
        verbose_name=_('Residues'),
        help_text=_('Amino acid residues that the respective codons code for')
    )
    start_codons = ArrayField(
        base_field=BooleanField(
            null=False,
            blank=False,
            unique=False,
            editable=True,
            db_index=False,
        ),
        size=64,
        null=False,
        blank=False,
        default=list,
        unique=False,
        editable=True,
        db_index=False,
        db_column='start_codons',
        verbose_name=_('Start Codons'),
        help_text=_('Whether or not the respective codon can act as a start codon')
    )
    stop_codons = ArrayField(
        base_field=BooleanField(
            null=False,
            blank=False,
            unique=False,
            editable=True,
            db_index=False,
        ),
        size=64,
        null=False,
        blank=False,
        default=list,
        unique=False,
        editable=True,
        db_index=False,
        db_column='stop_codons',
        verbose_name=_('Stop Codons'),
        help_text=_('Whether or not the respective codon can act as a stop codon')
    )
    # METHODS -----------------------------------------------------------------
    def __str__(self):
        """Return string representation"""
        return str(_(f'Translation Table: {self.genetic_code_id}'))

    def get_absolute_url(self):
        """Return URL"""
        return reverse_lazy('biology:genetic-code-detail', kwargs={'pk': self.pk},)

    def get_or_fetch(self):
        """Get object from database if it exists. Otherwise fetch from online"""
        try:
            entry = GeneticCode.objects.get(genetic_code_id=self.genetic_code_id)
        except GeneticCode.DoesNotExist:
            from twisted.internet import reactor, defer
            from scrapy.crawler import CrawlerRunner
            from scrapy.utils.project import get_project_settings
            from biocrawler.spiders import NCBIGeneticCodeSpider
            settings = get_project_settings()
            runner = CrawlerRunner(settings)
            @defer.inlineCallbacks
            def crawl():
                yield runner.crawl(NCBIGeneticCodeSpider)
                reactor.stop()
            crawl()
            reactor.run()
            entry = GeneticCode.objects.get(genetic_code_id=self.genetic_code_id)
        finally:
            return entry
