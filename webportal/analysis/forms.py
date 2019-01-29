from django import forms

class FastqForm(forms.Form):
    fastq1 = forms.FileField()
