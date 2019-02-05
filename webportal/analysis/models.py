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
    upload_name = models.CharField(max_length=200, blank=False)
    genome = models.CharField(max_length=200, choices=GENOME_CHOICES)
    organism = models.CharField(max_length=200)
    # conditions = models.CharField(max_length=50, blank=True, null=True)
    # no_replicates = models.PositiveSmallIntegerField()

    # def __str__(self):
    #     return self.upload_name


class Conditions(models.Model):
    session = models.ForeignKey(Session, on_delete=models.PROTECT)
    conditions = models.CharField(max_length=50, blank=False)
    no_replicates = models.PositiveSmallIntegerField(blank=False, default=1)

    def __str__(self):
        return self.conditions


class Samples(models.Model):
    LIBTYPE_CHOICES = (
        ("PE", "Paired_end"),
        ("SG", "Single")
    )

    session = models.ForeignKey(Session, on_delete=models.PROTECT)
    condition = models.ForeignKey(Conditions, on_delete=models.PROTECT)
    # replicate = models.PositiveSmallIntegerField()
    libtype = models.CharField(max_length=200, choices=LIBTYPE_CHOICES, blank=True, null=True)
    read_1 = models.FileField(upload_to='data/', blank=False)
    read_2 = models.FileField(upload_to='data/', blank=True, null=True)

    # def __str__(self):
    #     return self.sessions



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

# inspect available object methods
# dir(Product.objects.get(id=1))
