from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.
def index(request):
    return HttpResponse("Yo!")

def fileUpload(request):
    return HttpResponse("This is the page to upload files")
