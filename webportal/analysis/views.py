from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from .forms import SessionForm, WorkflowForm, SamplesForm, ConditionsForm
from analysis.models import Session, Workflow, Samples, Conditions
from . import models
from django.forms import modelformset_factory
# from django.views.generic import TemplateView
# from django.core.files.storage import FileSystemStorage
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)


###### Class based views ######
# Session
class cbv_view(View):
    def get(self, request):
            return HttpResponse('im not sure about these new fangled class views...')
            # return render(request,'index.html')

class SessionListView(ListView):
    # by default context object used to access parameters is lower(model)_list
    # context_object_name = 'session'
    model = models.Session

class SessionDetailView(DetailView):
    context_object_name = 'session_detail'
    model = models.Session
    queryset = Session.objects.all()
    # template_name = 'analysis/session_detail.html'
    # defining custom grouping field in urls.py we must then manually extract the key from self.kwargs
    def get_object(self):
        print(f'\nself.kwargs: {self.kwargs}')
        pk = self.kwargs.get('pk') # passes key 'pk' into kwargs dict to extract value
        print(f'pk: {pk}')
        return get_object_or_404(Session, pk=pk)


class SessionCreateView(CreateView):
    fields = ('organism','genome','fasta_file', 'annotation_file',)
    model = models.Session


class SessionUpdateView(UpdateView):
    fields = ('organism','genome','fasta_file', 'annotation_file',)
    model = models.Session


class SessionDeleteView(DeleteView):
    model = models.Session
    success_url = reverse_lazy("analysis:session_list")


# Conditions
class ConditionsListView(ListView):
    # context_object_name = 'conditions'
    model = models.Conditions

class ConditionsDetailView(DetailView):
    context_object_name = 'conditions_detail'
    model = models.Conditions
    template_name = 'analysis/conditions_detail.html'

# removing the session field from conditions create and using primary key defined by session
# in the url to pass to the conditions model. not working yet.
# NOT NULL constraint failed: analysis_conditions.session_id
class ConditionsCreateView(CreateView):
    template_name = 'analysis/conditions_form.html'
    # queryset = Conditions.objects.all()
    model = models.Conditions
    fields = ('session', 'conditions', 'no_replicates',)

    # def get(self, request, pk):
    #     form = ConditionsForm
    #     pk = self.kwargs.get('pk')
    #     print(f'\nPrimary Key: {pk}')
    #     return render(request, self.template_name, {'form':form})
    #
    # def post(self, request, pk):
    #     form = ConditionsForm(request.POST)
    #     if form.is_valid():
    #         print(form.cleaned_data['session'])
    #
    #         pk = self.kwargs.get('pk')
    #         post = form.save(commit=False)
    #         post.session = 'session' + str(pk)
    #         print(post.session)
    #         post.save()
    #         # print(f'\nRaw Post Obj: {post}')
    #         con = form.cleaned_data['conditions']
    #         no_rep = form.cleaned_data['no_replicates']
    #         # print(no_rep)
    #         return redirect('analysis:session_list')

# Cannot assign "'session1'": "Conditions.session" must be a "Session" instance.


class ConditionsUpdateView(UpdateView):
    fields = ('conditions','no_replicates',)
    model = models.Conditions

class ConditionsDeleteView(DeleteView):
    context_object_name = 'condition'
    model = models.Conditions
    success_url = reverse_lazy("analysis:session_list")


# Samples
class SamplesListView(ListView):
    # context_object_name = 'conditions'
    model = models.Samples

class SamplesDetailView(DetailView):
    context_object_name = 'samples_detail'
    model = models.Samples
    template_name = 'analysis/samples_detail.html'

class SamplesCreateView(CreateView):
    # fields = ('condition','libtype', 'read_1','read_2',)
    # model = models.Samples
    template_name = 'analysis/samples_form.html'
    form_class = SamplesForm
    queryset = Samples.objects.all()

    def form_valid(self, form):
        print(form.cleaned_data)
        return super().form_valid(form)

class SamplesUpdateView(UpdateView):
    fields = ('libtype','read_1','read_2',)
    model = models.Samples

class SamplesDeleteView(DeleteView):
    context_object_name = 'sample'
    model = models.Samples
    success_url = reverse_lazy("analysis:samples_list")



# Workflow
class WorkflowListView(ListView):
    # context_object_name = 'conditions'
    model = models.Workflow

class ConditionsDetailView(DetailView):
    context_object_name = 'conditions_detail'
    model = models.Conditions
    template_name = 'analysis/conditions_detail.html'

# removing the session field from conditions create and using primary key defined by session
# in the url to pass to the conditions model. not working yet.
# NOT NULL constraint failed: analysis_conditions.session_id
class ConditionsCreateView(CreateView):
    template_name = 'analysis/conditions_form.html'
    # queryset = Conditions.objects.all()
    model = models.Conditions
    fields = ('session', 'conditions', 'no_replicates',)

    # def get(self, request, pk):
    #     form = ConditionsForm
    #     pk = self.kwargs.get('pk')
    #     print(f'\nPrimary Key: {pk}')
    #     return render(request, self.template_name, {'form':form})
    #
    # def post(self, request, pk):
    #     form = ConditionsForm(request.POST)
    #     if form.is_valid():
    #         print(form.cleaned_data['session'])
    #
    #         pk = self.kwargs.get('pk')
    #         post = form.save(commit=False)
    #         post.session = 'session' + str(pk)
    #         print(post.session)
    #         post.save()
    #         # print(f'\nRaw Post Obj: {post}')
    #         con = form.cleaned_data['conditions']
    #         no_rep = form.cleaned_data['no_replicates']
    #         # print(no_rep)
    #         return redirect('analysis:session_list')

# Cannot assign "'session1'": "Conditions.session" must be a "Session" instance.


class ConditionsUpdateView(UpdateView):
    fields = ('conditions','no_replicates',)
    model = models.Conditions

class ConditionsDeleteView(DeleteView):
    context_object_name = 'condition'
    model = models.Conditions
    success_url = reverse_lazy("analysis:session_list")












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




def detail(request, venue_id):
    session = get_object_or_404(Venue, pk=venue_id)
    return render(request, 'analysis/conditions_list.html', {'session': session})
