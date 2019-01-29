from django.contrib import admin

# Register your models here.
from .models import Session, Workflow, Samples, Post #import models

admin.site.register(Session)
admin.site.register(Workflow)
admin.site.register(Samples)
admin.site.register(Post)

# username: teamrnaseq
#pass: passwordabc # note: NOT for production!
