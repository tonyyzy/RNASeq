from django import forms
from django.forms import formset_factory


class UserCreationForm(forms.Form):
    username = forms.CharField(max_length=30)
    password1 = forms.CharField(widget=forms.PasswordInput)
    password2 = forms.CharField(widget=forms.PasswordInput)

    def clean_password2(self):
        password1 = self.cleaned_data.get("password1", "")
        password2 = self.cleaned_data["password2"]
        if password1 != password2:
            raise forms.ValidationError('passwords did not match')
        return password2


class BookForm(forms.Form):
    name = forms.CharField(
            max_length=30,
            label='Book Name',
            widget=forms.TextInput(attrs={
            'class': 'form-control',
            'placeholder': 'Enter Book Name here'
            })
        )
    author = forms.CharField(
            max_length=30,
            label='author Name',
            widget=forms.TextInput(attrs={
            'class': 'form-control',
            'placeholder': 'Enter author Name here'
            })
        )

BookFormset = formset_factory(BookForm, extra=2)
#
# ,
# author = forms.CharField(
#     label='Author Name',
#     widget=forms.TextInput(attrs={
#         'class': 'form-control',
#         'placeholder': 'Enter Author here'
#     })
