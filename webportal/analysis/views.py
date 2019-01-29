from django.shortcuts import render
from django.http import HttpResponse
from .forms import *


# Create your views here.
# def index(request):
#     return HttpResponse("Yo!")


def index(request):
    return HttpResponse("you should not be reading this. it means something broke")
    # my_dict = {'insert': 'inserted from analysis'}
    # return render(request, 'templates/index.html', context=my_dict)

def test(request):
    return HttpResponse("testing page for debug")
    # my_dict = {'insert': 'inserted from analysis'}
    # return render(request, 'templates/test.html', context=my_dict)


def get_files(request):
    # return HttpResponse("post me post me post me post me post me post me")
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = FileSubmission(request.POST, request.FILES) # check whether it's valid:
        if form.is_valid():
            text = form.cleaned_data['post'] # assign clean data to text variable

        args = {'form': form, 'text': text}
        return render(request, 'file_submit.html', args)

    # if a GET (or any other method) we'll create a blank form
    else:
        form = FileSubmission()
    return render(request, 'file_submit.html', {'form': form})





def session_form_upload(request):
    if request.method == 'POST':
        form = SessionForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            return redirect('home')
    else:
        form = SessionForm()
    return render(request, 'session.html', {
        'form': form
    })
