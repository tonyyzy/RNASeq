from django.shortcuts import render, redirect
from django.http import HttpResponse
from analysis.models import Product, Fastq
from .form import ProductForm, FastqForm
from django.forms import modelformset_factory

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

def thanks_view(request):
    # return HttpResponse("pages/views homes_view return.")
    return render(request, 'thanks.html', {})

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


# def create_view(request):
#     if request.method == 'POST':
#         form = FastqForm(request.POST, request.FILES)
#         if form.is_valid():
#             form.save()
#     else:
#         form = FastqForm()
#     args = {'form': form}
#     return render(request, 'create.html', args)

# create with form factory
def create_view(request):
    fastq_formset = modelformset_factory(Fastq, fields=('name', 'fastq', 'library', 'condition'), extra=2)

    if request.method == 'POST':
        form = fastq_formset(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            return redirect('/thanks')
    else:
        form = fastq_formset(queryset=Product.objects.none())
    args = {'form': form}
    return render(request, 'create.html', args)



# def create_factory(request):
#     formSet = modelformset_factory(Product, fields=('title', 'description', 'price', 'summary', 'featured'), extra=4)
#
#     if request.method == 'POST':
#
#         form = formSet(request.POST)
#
#         instances = form.save(commit=False)
#         for instance in instances:
#             instance.save()
#
#     form = formSet(queryset=Product.objects.none()) # queryset used to prevent current db display
#     args = {'form': form}
#     return render(request, 'create.html', args)
