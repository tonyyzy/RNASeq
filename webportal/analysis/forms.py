from django import forms
from django.forms import ModelForm
from analysis.models import Session, Workflow, Samples


class SessionForm(forms.ModelForm):
    class Meta:
        model = Session
        # fields = ['GENOME_CHOICES', 'genome', 'organism', 'status', 'no_conditions', 'no_replicates']
        fields='__all__'

class SamplesForm(forms.ModelForm):
    class Meta:
        model = Samples
        # fields = ['GENOME_CHOICES', 'genome', 'organism', 'status', 'no_conditions', 'no_replicates']
        fields='__all__'


class WorkflowForm(forms.ModelForm):
    class Meta:
        model = Workflow
        # fields = ['GENOME_CHOICES', 'genome', 'organism', 'status', 'no_conditions', 'no_replicates']
        fields='__all__'
