from django import forms
from django.forms import ModelForm
from analysis.models import Session, Workflow, Samples, Conditions, The_Debug, Genome
from django.core.exceptions import ValidationError

class SessionSearchForm(forms.Form):
    user_session = forms.CharField(max_length=100)


class GenomeForm(forms.ModelForm):
    class Meta:
        model = Genome
        fields='__all__'


class SessionForm(forms.ModelForm):

    class Meta:
        model = Session
        fields = ['genome_index','select_genome','organism', 'salmon', 'fasta_dna_file', 'fasta_cdna_file', 'gtf_file', 'status']
        # fields='__all__'
        widgets={
            'genome_index': forms.Select(attrs={
                'class':'form-control',
                }),
            'select_genome': forms.Select(attrs={
                'class':'form-control',
                }),
            'organism': forms.TextInput(attrs={
                'class':'form-control',
                'placeholder': 'enter organism name here...'
                }),
            'salmon': forms.NullBooleanSelect(attrs={
                'class':'form-control',
                }),
}



class SessionSubmitForm(forms.ModelForm):

    class Meta:
        model = Session
        fields = []


class ConditionsForm(forms.ModelForm):
    class Meta:
        model = Conditions
        fields = ['conditions', 'no_replicates']
        # fields='__all__'
        widgets={
            'conditions': forms.TextInput(attrs={
                'class':'form-control',
                'placeholder': 'enter condition here...'
                }),
        }



class SamplesForm(forms.ModelForm):
    class Meta:
        model = Samples
        # fields = ['condition', 'libtype', 'read_1']
        fields = ['libtype', 'read_1', 'read_2', 'accession']
        # fields = ['session', 'condition', 'libtype', 'read_1', 'read_2', 'accession']
        # fields='__all__'


class WorkflowForm(forms.ModelForm):
    class Meta:
        model = Workflow
        fields = ['mapper', 'assembler', 'analysis', 'status']
        # fields='__all__'



class DebugForm(forms.ModelForm):
    debug_title = forms.CharField(
            max_length=100,
            widget=forms.TextInput(
                attrs={'class':'form-control', 'placeholder':'debug_title'}
            )
        )

    debug_body = forms.CharField(
            max_length=100,
            widget=forms.Textarea(
                attrs={'class':'form-control', 'placeholder':'debug_body'}
            )
        )

    class Meta:
        model = The_Debug
        fields = ['field_one', 'field_two', 'field_three']
        # fields='__all__'
        widgets = {
            'field_one': forms.TextInput(
                attrs={
                    'class':'form-control',
                    'placeholder':'field_one_input',
                    }),
            'field_two': forms.TextInput(
                attrs={
                    'class':'form-control',
                    'placeholder':'field_two_input',
                    }),
            'field_three': forms. FileInput(
                attrs={
                    'placeholder':'field_three_input',
                    })
        }
