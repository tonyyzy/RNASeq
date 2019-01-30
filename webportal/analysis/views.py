from django.shortcuts import render
from django.http import HttpResponse
from analysis.forms import homeForm
from .forms import FileSubmission
from .models import Product
# from django.views.generic import TemplateView



def my_view(request):
    # form = Product()
    obj = Product.objects.get(id=1)
    # print(posts)
    args = {
        'title': obj.title,
        'description': obj.description,
        'price': obj.price}
    # args = {'object': obj}
    return render(request, 'output/detail.html', args)
    # return render(request, {'form': form})


# class analysisView(TemplateView):
#     template_name = 'analysis/templates/analysis'
#
#     def get(self, request):
#         form = homeForm()
#         return render(request, self.template_name, {'form': form})


def index(request):
    # return HttpResponse("you should not be reading this. it means something broke")
    my_dict = {'insert': 'content rendered from view.py'}
    return render(request, 'analysis/index.html', context=my_dict) # renders page from templates dir

def test(request):
    return HttpResponse("testing page for debug")
    # my_dict = {'insert': 'inserted from analysis'}
    # return render(request, 'templates/test.html', context=my_dict)


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

#
# def session_form_upload(request):
#     if request.method == 'POST':
#         form = SessionForm(request.POST, request.FILES)
#         if form.is_valid():
#             form.save()
#             return redirect('home')
#     else:
#         form = SessionForm()
#     return render(request, 'session.html', {
#         'form': form
#     })
