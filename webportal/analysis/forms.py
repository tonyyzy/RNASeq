from django import forms
from django.forms import ModelForm
from analysis.models import Session, Workflow, Samples, Conditions
from django.core.exceptions import ValidationError

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
        # fields = ['condition', 'libtype', 'read_1']
        fields = ['libtype', 'read_1', 'read_2', 'accession']
        # fields = ['session', 'condition', 'libtype', 'read_1', 'read_2', 'accession']
        # fields='__all__'

# attempted django client side form validation not working currently.
        # def clean(self):
            # data = super(SamplesForm, self).clean()
            # print(f'\n{data}\n')
        #
        # def clean(self):
        #     lib = self.cleaned_data.get('libtype')
        #     if lib == 'PE':
        #         msg = forms.ValidationError("This field is required.")
        #         self.add_error('shipping_destination', msg)
        #     else:
        #     # Keep the database consistent. The user may have
        #     # submitted a shipping_destination even if shipping
        #     # was not selected
        #     # self.cleaned_data['shipping_destination'] = ''
        #         pass
        #     return self.cleaned_data
                # raise ValidationError('Telephone is required')
            # # read_2 = self.cleaned_data.get('read_2')
            # if not libtype == 'PE':
            #     raise forms.ValidationError('invalid!')
            # return libtype


class WorkflowForm(forms.ModelForm):
    class Meta:
        model = Workflow
        fields = ['index', 'mapper', 'assembler', 'analysis', 'status']
        # fields='__all__'



# pip install --upgrade django-crispy-forms
# def clean(self):
#     shipping = self.cleaned_data.get('shipping')
#
#     if shipping:
#         msg = forms.ValidationError("This field is required.")
#         self.add_error('shipping_destination', msg)
#     else:
#         # Keep the database consistent. The user may have
#         # submitted a shipping_destination even if shipping
#         # was not selected
#         self.cleaned_data['shipping_destination'] = ''
#
#     return self.cleaned_data
