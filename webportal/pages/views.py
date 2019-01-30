from django.shortcuts import render
from django.http import HttpResponse
from analysis.models import Product

# Create your views here.
def home_view(request):
    # return HttpResponse("pages/views homes_view return.")
    return render(request, 'index.html', {})

def test_view(request):
    return HttpResponse("pages/views.py test_view return")

def about_view(request):
    # return HttpResponse("pages/views homes_view return.")
    return render(request, 'about.html', {})

def analysis_view(request):
    # return HttpResponse("pages/views homes_view return.")
    return render(request, 'analysis.html', {})

def detail_view(request):
    # return HttpResponse("pages/views homes_view return.")
    return render(request, 'detail.html', {})

def db_view(request):
    obj = Product.objects.get(id=1)
    args = {'object': obj}
    return render(request, 'output/detail.html', args)
    
    # for loop version
    # args = {
    #     'title': obj.title,
    #     'description': obj.description,
    #     'price': obj.price,
    #     'summary': obj.summary,
    #     'price': obj.featured,
    #     }
    # return render(request, 'output/detail.html', {'args': args})



def upload_view(request):
    # return HttpResponse("pages/views homes_view return.")
    return render(request, 'detail.html', {})
