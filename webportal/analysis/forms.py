from django import forms
from django.forms import ModelForm
from analysis.models import Session, Workflow, Samples, Conditions


class SessionSearchForm(forms.Form):
    user_session = forms.CharField(max_length=100)

class SessionForm(forms.ModelForm):
    class Meta:
        model = Session
        fields = ['organism','genome','fasta_file', 'annotation_file']
        # fields='__all__'

class ConditionsForm(forms.ModelForm):
    class Meta:
        model = Conditions
        fields = ['conditions', 'no_replicates']
        # fields = ['session', 'conditions', 'no_replicates']
        # fields='__all__'


class SamplesForm(forms.ModelForm):
    class Meta:
        model = Samples
        fields = ['condition', 'libtype', 'read_1', 'read_2']
        # fields='__all__'


class WorkflowForm(forms.ModelForm):
    class Meta:
        model = Workflow
        fields = ['index', 'mapper', 'assembler', 'analysis', 'status']
        # fields='__all__'



# pip install --upgrade django-crispy-forms
