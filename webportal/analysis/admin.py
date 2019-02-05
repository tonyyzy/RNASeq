from django.contrib import admin

# Register your models here.
from .models import Session, Workflow, Samples #import models

admin.site.register(Session)
admin.site.register(Workflow)
admin.site.register(Samples)

# username: teamrnaseq
#pass: passwordabc # note: NOT for production!
