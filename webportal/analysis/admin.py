from django.contrib import admin

# Register your models here.
from .models import Session, Workflow, Samples, Product, Fastq #import models

admin.site.register(Session)
admin.site.register(Workflow)
admin.site.register(Samples)
admin.site.register(Product)
admin.site.register(Fastq)

# username: teamrnaseq
#pass: passwordabc # note: NOT for production!
