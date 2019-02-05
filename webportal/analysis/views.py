from django.shortcuts import render, redirect
from django.http import HttpResponse
from .forms import SessionForm, WorkflowForm, SamplesForm, ConditionsForm
from analysis.models import Session, Workflow, Samples
from django.forms import modelformset_factory
# from django.views.generic import TemplateView
from django.core.files.storage import FileSystemStorage



def home_view(request):
    # return HttpResponse("what hath god wrought")
    return render(request, 'analysis/home.html', {})


def session_view(request):
    # return HttpResponse("what hath god wrought")
    if request.method == 'POST':
        form = SessionForm(request.POST, request.FILES)
        if form.is_valid():
            print('bound session form posted \n')
            form.save()
            # return render(request, 'analysis/upload_samples.html', {})
            return redirect('analysis:upload_session')
    else:
        form = SessionForm()
    return render(request, 'analysis/upload_session.html', {'form': form})


def conditions_view(request):
    # return HttpResponse("what hath god wrought")
    if request.method == 'POST':
        form = ConditionsForm(request.POST, request.FILES)
        if form.is_valid():
            print('bound conditions form posted \n')
            form.save()
            # return render(request, 'analysis/upload_samples.html', {})
            return redirect('analysis:upload_conditions')
    else:
        form = ConditionsForm()
    return render(request, 'analysis/upload_conditions.html', {'form': form})


def samples_view(request):
    # sample_formset = modelformset_factory(Samples, fields=('__all__'), extra=1)
    if request.method == 'POST':
        # form = sample_formset(request.POST, request.FILES)
        form = SamplesForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            print('bound samples form posted \n')
            # context = {'form': form}
            # return render(request, 'analysis/upload_workflow.html', context)
            return redirect('analysis:upload_samples')
    else:
        pass
        # form = sample_formset()
        # form = sample_formset(queryset=Samples.objects.filter(session=7))
        # form = sample_formset(queryset=Samples.objects.none())
        form = SamplesForm()
    return render(request, 'analysis/upload_samples.html', {'form': form})



def workflow_view(request):
    if request.method == 'POST':
        form = WorkflowForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            return render(request, 'thanks.html', {})
    else:
        form = WorkflowForm()
    return render(request, 'analysis/upload_workflow.html', {'form': form})


def samples_list_view(request):
    # return HttpResponse("post me post me post me post me post me post me")
    samples = Samples.objects.all()
    context = {'samples': samples}
    return render(request, 'analysis/samples_list.html', context)




# FILE SUBMISSION EXAMPLE
def get_files(request):
    # return HttpResponse("post me post me post me post me post me post me")
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = FileSubmission(request.POST, request.FILES) # check whether it's valid:
        # form = ModelFormWithFileField(request.POST, request.FILES)
        if form.is_valid():
            user = request.user
            post = form.save(commit=False)
            post.user = request.user
            post.save()
            text = form.cleaned_data['post'] # assign clean data to text variable

        args = {'form': form, 'text': text, 'user': user}
        return render(request, 'analysis/file_submit.html', args)

    # if a GET (or any other method) we'll create a blank form
    else:
        form = FileSubmission()
    return render(request, 'analysis/file_submit.html', {'form': form})
