from django.shortcuts import render
from django.http import HttpResponse
from analysis.forms import homeForm
from .forms import SessionForm, SamplesForm, WorkflowForm
# from django.views.generic import TemplateView



def test_view(request):
    # return HttpResponse("what hath god wrought")
    my_dict = {'insert': 'content rendered from analysis/views.py'}
    return render(request, 'analysis/index.html', context=my_dict) # renders page from templates dir


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


# def session_view(request):
#     if request.method == 'POST':
#         form = SessionForm(request.POST)
#         if form.is_valid():
#             form.save()
#             return render(request, 'thanks.html', {})
#     else:
#         form = SessionForm()
#     return render(request, 'analysis/upload.html', {'form': form})
#
#

def samples_view(request):
    if request.method == 'POST':
        form = SamplesForm(request.POST, request.FILES)
        if form.is_valid():
            form.save('/data/fastq')
            return render(request, 'thanks.html', {})
    else:
        form = SamplesForm()
    return render(request, 'analysis/upload.html', {'form': form})

# def samples_view(request):
#     if request.method == 'POST':
#         form = WorkflowForm(request.POST)
#         if form.is_valid():
#             form.save()
#             return render(request, 'thanks.html', {})
#     else:
#         form = WorkflowForm()
#     return render(request, 'analysis/upload.html', {'form': form})



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
