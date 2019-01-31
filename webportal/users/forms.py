from django import forms
from django.contrib.auth.models import User
from django.contrib.auth.forms import UserCreationForm


class userRegistrationForm(UserCreationForm):
    email = forms.EmailField()

    class Meta:
        model = User
        fields = [
        'first_name',
        'last_name',
        'username',
        'email',
        'password1',
        'password2',
        ]
