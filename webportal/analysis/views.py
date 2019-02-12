from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from .forms import SessionSearchForm, SessionForm, WorkflowForm, SamplesForm, ConditionsForm
from analysis.models import Session, Workflow, Samples, Conditions
from . import models
from django.forms import modelformset_factory
from django.db.models import Q
# from django.views.generic import TemplateView
# from django.core.files.storage import FileSystemStorage
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)


###### Class based views ######
# Session
class SessionIndexView(View):
    def get(self, request):
        # return HttpResponse('Session Index View')
        sessions = Session.objects.all()
        # form = SessionSearchForm
        # context = {'form': form}
        query = self.request.GET.get('q')
        if query:
            query = query.replace('-','')
            print(f'cleaned query: {query}')
            query_set_list = sessions.filter(
                            Q(identifier__icontains=query)|
                            Q(organism__icontains=query)
                            ).distinct()
            # query_set_list = query_set_list.first()

            print(query_set_list)
            context = {'query_set_list':query_set_list}
            return render(request, 'analysis/index.html', context)
        return render(request, 'analysis/index.html')

    def post(self, request):
        form = SessionSearchForm(request.POST)
        # context = {'form': form, 'sessions': sessions}
        if form.is_valid():
            user_session = form.cleaned_data['user_session']
            # print(user_session)
            context = {'form': form, 'sessions': sessions, 'user_session':user_session}
            return render(request, 'analysis/index.html', context)
            # return redirect('analysis:session_index')
        context = {'form': form, 'sessions': sessions}
        return render(request, 'analysis/index.html', context)

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


class ConditionsCreateView(CreateView):
    template_name = 'analysis/conditions_form.html'
    # queryset = Conditions.objects.all()

    def get(self, request, pk):
        form = ConditionsForm
        context = {'form':form}
        return render(request, self.template_name, context)
        print('\nFORM WAS VALID: POST.SESSION:')

    def post(self, request, pk):
        form = ConditionsForm(request.POST)
        if form.is_valid():
            post = form.save(commit=False)
            sessions = Session.objects.all()
            adjusted_pk = self.kwargs.get('pk')-1 # python index 1st session from 0
            # print(f'\n{sessions[adjusted_pk]}') # returns a session object NOT a string!
            post.session = sessions[adjusted_pk]
            post.save()
            return redirect('analysis:conditions_detail', pk)
        return render(request, self.template_name, {'form':form})



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

class WorkflowDetailView(DetailView):
    context_object_name = 'workflow_detail'
    model = models.Workflow
    template_name = 'analysis/workflow_detail.html'

class WorkflowCreateView(CreateView):
    template_name = 'analysis/workflow_form.html'
    # queryset = Conditions.objects.all()
    model = models.Workflow
    fields = ('session', 'index', 'mapper', 'assembler', 'analysis', 'status',)

class WorkflowUpdateView(UpdateView):
    fields = ('mapper', 'assembler', 'analysis', 'status',)
    model = models.Workflow

class WorkflowDeleteView(DeleteView):
    context_object_name = 'workflow'
    model = models.Workflow
    success_url = reverse_lazy("analysis:workflow_list")












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
