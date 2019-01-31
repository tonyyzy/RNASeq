from django import forms
from analysis.models import Product, Fastq


class ProductForm(forms.ModelForm):
    class Meta:
        model = Product
        fields = [
            'title',
            'description',
            'price',
            'summary',
            'featured',
        ]
class FastqForm(forms.ModelForm):
    class Meta:
        model = Fastq
        fields = [
            # 'identifier' #not required as non editable field
            'name',
            'fastq',
            'library',
            'condition',
            # 'created_at',
        ]
