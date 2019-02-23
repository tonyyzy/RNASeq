from django.db import models
import uuid
from django.contrib.auth.models import User
from django.conf import settings
from django.urls import reverse


class Session(models.Model):
    identifier = models.UUIDField(default=uuid.uuid4, editable=False, unique=True)
    GENOME_CHOICES = (
        ("PreIndex", "Preindexed_Genome"),
        ("OwnGFF", "User_Provided_Annotation"),
        ("OwnFASTA", "User_Provided_Genome"),
    )
    organism = models.CharField(max_length=200)
    genome = models.CharField(max_length=200, choices=GENOME_CHOICES)
    fasta_file = models.FileField(upload_to='data/', blank=True, null=True)
    annotation_file = models.FileField(upload_to='data/', blank=True, null=True)

    def get_absolute_url(self): # provides a default if Session is called from views.py without a specified reverse or redirect
        return reverse('analysis:session_detail', kwargs={'session_slug':self.identifier})

    def __str__(self):
        return 'session' + str(self.pk)


class Conditions(models.Model):
    session = models.ForeignKey(Session, on_delete=models.PROTECT, related_name='conditions_fk')
    conditions = models.CharField(max_length=50, blank=False)
    no_replicates = models.PositiveSmallIntegerField(blank=False, default=1)

    def get_absolute_url(self):
        return reverse('analysis:session_detail', kwargs={'pk':self.pk})

    def __str__(self):
        return self.conditions


class Samples(models.Model):
    LIBTYPE_CHOICES = (
        ("PE", "Paired_end"),
        ("SG", "Single")
    )
    session = models.ForeignKey(Session, on_delete=models.PROTECT, related_name='samples_fk')
    condition = models.ForeignKey(Conditions, on_delete=models.PROTECT)
    libtype = models.CharField(max_length=200, choices=LIBTYPE_CHOICES, blank=False, null=False)
    read_1 = models.FileField(upload_to='data/', blank=False, null=False)
    read_2 = models.FileField(upload_to='data/', blank=True, null=True)
    accession = models.CharField(max_length=200, blank=False, null=False)

    def get_absolute_url(self):
        return reverse('analysis:session_detail', kwargs={'slug':self.slug})


class Workflow(models.Model):
    INDEX_CHOICES = (
        ("STAR", "STARAligner"),
        ("HISAT2", "HISAT2"),
    )
    MAPPER_CHOICES = (
        ("STAR", "STARAligner"),
        ("HISAT2", "HISAT2"),
    )
    ASSEMLBER_CHOICES = (
        ("STRINGTIE", "STRINGTIE"),
    )
    ANALYSIS_CHOICES = (
        ('DESEQ2', 'DESEQ2'),
        ('DEXEQ', 'DEXEQ'),
        ('HISAT2', 'HISAT2'),
        ('HTSEQ', 'HTSEQ'),
        ('PREPDE', 'PREPDE'),
        ('SAMTOOLS', 'SAMTOOLS'),
        ('STAR', 'STAR'),
        ('STRINGTIE', 'STRINGTIE'),
        ('MISO', 'MISO'),
        ('SALMON', 'SALMON'),
        ('DESEQ', 'DESEQ'),
        ('CUFFLINKS', 'CUFFLINKS'),
    )
    session = models.ForeignKey(Session, on_delete=models.PROTECT, related_name='workflow')
    index = models.CharField(max_length=200, choices=INDEX_CHOICES)
    mapper = models.CharField(max_length=200, choices=MAPPER_CHOICES)
    assembler = models.CharField(max_length=200, choices=ASSEMLBER_CHOICES, blank=True)
    analysis = models.CharField(max_length=200, choices=ANALYSIS_CHOICES)
    status = models.BooleanField(default=False, null=False)

    def get_absolute_url(self):
        return reverse('analysis:session_detail', kwargs={'pk':self.pk})
