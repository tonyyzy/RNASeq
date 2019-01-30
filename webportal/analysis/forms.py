from django import forms
from analysis.models import Session, Workflow, Samples, Product


class FileSubmission(forms.ModelForm):
    product = forms.CharField()
    class Meta:
        model = Product
        fields = ('product',)

class homeForm(forms.Form):
    product = forms.CharField()
