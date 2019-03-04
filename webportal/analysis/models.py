from django.db import models
import uuid
from django.contrib.auth.models import User
from django.conf import settings
from django.urls import reverse
from django.core.files.storage import FileSystemStorage
import os

data_root = FileSystemStorage(location=settings.DATA_DIR)
genome_index_root = FileSystemStorage(location=settings.GENOME_INDEX_DIR)


class Genome(models.Model):
    organism = models.CharField(max_length=500)
    source = models.CharField(max_length=500)
    version = models.CharField(max_length=500)
    fasta_dna_file = models.CharField(max_length=500, blank=False, null=False)
    fasta_cdna_file = models.CharField(max_length=500, blank=False, null=False)
    gtf_file = models.CharField(max_length=500, blank=False, null=False)
    star = models.CharField(max_length=500)
    salmon = models.CharField(max_length=500)
    hisat2 = models.CharField(max_length=500)

    def __str__(self):
        return self.organism


class Session(models.Model):

    def get_genome_path(self, filename):
        return os.path.join(str(self.identifier), 'genome', filename) # version with dash
        # return os.path.join(self.identifier.hex, filename) # version without dash

    GENOME_CHOICES = (
        ("pre_index", "Preindexed Genome"),
        ("user_provided", "Provide Own Index"),
    )

    identifier = models.UUIDField(default=uuid.uuid4, editable=False, unique=True)
    genome_index = models.CharField(max_length=200, choices=GENOME_CHOICES)
    select_genome = models.ForeignKey(Genome, on_delete=models.PROTECT, related_name='genome_fk', blank=True, null=True)
    organism = models.CharField(max_length=200, blank=True, null=True)
    salmon = models.BooleanField(blank=True, null=True)
    fasta_dna_file = models.FileField(storage=data_root, upload_to=get_genome_path, blank=True, null=True)
    fasta_cdna_file = models.FileField(storage=data_root, upload_to=get_genome_path, blank=True, null=True)
    gtf_file = models.FileField(storage=data_root, upload_to=get_genome_path, blank=True, null=True)
    status = models.BooleanField(default=False, blank=True, null=True)

    def get_absolute_url(self): # provides a default if Session is called from views.py without a specified reverse or redirect
        return reverse('analysis:session_detail', kwargs={'session_slug':self.identifier})

    def __str__(self):
        return 'session' + str(self.pk)


class Condition(models.Model):
    session = models.ForeignKey(Session, on_delete=models.PROTECT, related_name='conditions_fk')
    condition = models.CharField(max_length=50, blank=False)
    no_replicates = models.PositiveSmallIntegerField(blank=False, default=1)

    def get_absolute_url(self):
        return reverse('analysis:session_detail', kwargs={'pk':self.pk})

    def __str__(self):
        return self.condition


class Samples(models.Model):
    LIBTYPE_CHOICES = (
        ("PE", "Paired_end"),
        ("SG", "Single")
    )

    def get_fastq_path(self, filename):
        return os.path.join(str(self.session.identifier), 'fastq', filename) # version with dash
        # return os.path.join(self.identifier.hex, filename) # version without dash

    session = models.ForeignKey(Session, on_delete=models.PROTECT, related_name='samples_fk')
    condition = models.ForeignKey(Condition, on_delete=models.PROTECT)
    libtype = models.CharField(max_length=200, choices=LIBTYPE_CHOICES, blank=False, null=False)
    read_1 = models.FileField(storage=data_root, upload_to=get_fastq_path, blank=False, null=False)
    read_2 = models.FileField(storage=data_root, upload_to=get_fastq_path, blank=True, null=True)
    accession = models.CharField(max_length=200, blank=False, null=False)

    def get_absolute_url(self):
        return reverse('analysis:session_detail', kwargs={'slug':self.slug})


class Workflow(models.Model):
    # INDEX_CHOICES = (
    #     ("star", "STARAligner"),
    #     ("hisat2", "HISAT2"),
    #     ('salmon', 'SALMON'),
    # )
    MAPPER_CHOICES = (
        ("star", "STARAligner"),
        ("hisat2", "HISAT2"),
        ('salmonquant', 'SALMON'),
    )
    ASSEMLBER_CHOICES = (
        ("stringtie", "STRINGTIE"),
        ('cufflinks', 'CUFFLINKS'),
        ('miso', 'MISO'),
        ('htseq', 'HTSEQ'),
        ('featurecounts', 'FEATURECOUNTS'),
        ('salmoncount', 'SALMON')
    )
    ANALYSIS_CHOICES = (
        ('deseq2', 'DESEQ2'),
        ('dexseq', 'DEXSEQ'),
        ('htseq', 'HTSEQ'),
        ('miso', 'MISO'),
        ('cuffdiff', "CUFFDIFF"),
        ('edger', "EDGER"),
        ('ballgown', "BALLGOWN")
    )
    session = models.ForeignKey(Session, on_delete=models.PROTECT, related_name='workflow')
    # index = models.CharField(max_length=200, choices=INDEX_CHOICES)
    mapper = models.CharField(max_length=200, choices=MAPPER_CHOICES)
    assembler = models.CharField(max_length=200, choices=ASSEMLBER_CHOICES, blank=True)
    analysis = models.CharField(max_length=200, choices=ANALYSIS_CHOICES)
    status = models.BooleanField(default=False, null=False)

    def get_absolute_url(self):
        return reverse('analysis:session_detail', kwargs={'pk':self.pk})



class The_Debug(models.Model):

    def get_upload_path(self, filename):
        return os.path.join(str(self.identifier), filename) # version with dash
        # return os.path.join(self.identifier.hex, filename)

    FIELD_THREE_CHOICES = (
        ("choice1", "choice1"),
        ("choice2", "choice2"),
        ("choice3", "choice3"),
    )
    identifier = models.UUIDField(default=uuid.uuid4, editable=False, unique=True)
    field_one = models.CharField(max_length=200)
    # field_two = models.CharField(max_length=200)
    field_two = models.FileField(storage=data_root, upload_to=get_upload_path)
