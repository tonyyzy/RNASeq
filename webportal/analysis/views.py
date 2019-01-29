from django.shortcuts import render
from django.http import HttpResponse
from .forms import *

# Create your views here.
def index(request):
    return HttpResponse("Yo!")

def fileUpload(request):
    return HttpResponse("This is the page to upload files")

def get_fastq(request):
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = FastqForm(request.POST, request.FILES)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            return HttpResponse('Thank you!')

    # if a GET (or any other method) we'll create a blank form
    else:
        form = FastqForm()

    return render(request, 'fastq.html', {'form': form})
