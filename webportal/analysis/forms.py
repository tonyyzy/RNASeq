from django import forms
from analysis.models import Session, Workflow, Samples, Post


class FileSubmission(forms.ModelForm):
    post = forms.CharField()
    class Meta:
        model = Post
        fields = ('post',)



class SessionForm(forms.ModelForm):
    class Meta:
        model = Session
        fields = ('genome', 'organism')
