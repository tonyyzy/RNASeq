from django.db import models
import uuid
from django.contrib.auth.models import User
from django.conf import settings


class Session(models.Model):
    identifier = models.UUIDField(default=uuid.uuid4, editable=False, unique=True)
    GENOME_CHOICES = (
        ("PreIndex", "Preindexed_Genome"),
        ("OwnGFF", "User_Provided_Annotation"),
        ("OwnFASTA", "User_Provided_Genome"),
    )
    genome = models.CharField(max_length=200, choices=GENOME_CHOICES)
    organism = models.CharField(max_length=200)
    status = models.BooleanField(default=False, null=False)
    no_conditions = models.PositiveSmallIntegerField()
    no_replicates = models.PositiveSmallIntegerField()

    def __str__(self):
        return self.organism

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
    session = models.ForeignKey(Session, on_delete=models.PROTECT)
    index = models.CharField(max_length=200, choices=INDEX_CHOICES)
    mapper = models.CharField(max_length=200, choices=MAPPER_CHOICES)
    assembler = models.CharField(max_length=200, choices=ASSEMLBER_CHOICES, blank=True)
    analysis = models.CharField(max_length=200)
    status = models.BooleanField(default=False, null=False)

class Samples(models.Model):
    LIBTYPE_CHOICES = (
        ("PE", "Paired_end"),
        ("SG", "Single")
    )
    READ_CHOICES = (
        (1, "First"),
        (2, "Second")
    )
    session = models.ForeignKey(Session, on_delete=models.PROTECT)
    condition = models.CharField(max_length=200)
    # replicate = models.PositiveSmallIntegerField()
    libtype = models.CharField(max_length=200, choices=LIBTYPE_CHOICES, blank=True, null=True)
    read_1 = models.FileField(upload_to='data/', blank=False)
    # strand_1 = models.PositiveSmallIntegerField(choices=READ_CHOICES, default=1)
    read_2 = models.FileField(upload_to='data/', blank=True, null=True)
    # strand = models.PositiveSmallIntegerField(choices=READ_CHOICES, blank=True)

    # def __str__(self):
    #     return self.sessions

# inspect available object methods
# dir(Product.objects.get(id=1))
