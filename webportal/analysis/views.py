from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from .forms import SessionSearchForm, SessionForm, WorkflowForm, SamplesForm, ConditionsForm
from analysis.models import Session, Workflow, Samples, Conditions
from . import models
from django.db.models import Q
from django.contrib import messages
# from django.core.files.storage import FileSystemStorage
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)


# Session
class SessionIndexView(View):
    def get(self, request):
        sessions = Session.objects.all()
        query = self.request.GET.get('q')
        if query:
            print(f'\n{query}\n')
            # query_set = sessions.filter(Q(identifier__exact=query)
            try:
                session_query = Session.objects.get(identifier__exact=query)
                print(session_query)
                messages.success(request, 'success')
                context = {'session_query':session_query}
                return render(request, 'analysis/index.html', context)
            except:
                messages.warning(request, 'Session ID not found')
                context = {'invalid_query': query}
                return render(request, 'analysis/index.html', context)
        return render(request, 'analysis/index.html')



class SessionListView(ListView):
    # context_object_name = 'session'
    # by default context object used to access database object within template is is lower(model_name)_list
    model = models.Session


class SessionDetailView(DetailView):
    context_object_name = 'session_detail'
    # template_name = 'analysis/session_detail.html'
    model = models.Session

    def get_object(self):
        slug = self.kwargs.get('session_slug')
        return get_object_or_404(Session, identifier=slug)


class SessionCreateView(CreateView): # CreateView expects template lower(model_name)_form.html
    fields = ('organism','genome','fasta_file', 'annotation_file',)
    model = models.Session


class SessionUpdateView(UpdateView):
    fields = ('organism','genome','fasta_file', 'annotation_file',)
    model = models.Session


class SessionDeleteView(DeleteView):
    model = models.Session
    # success_url = reverse_lazy("analysis:session_list")


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

    def get(self, request, session_slug):
        form = ConditionsForm
        context = {'form':form}
        return render(request, self.template_name, context)

    def post(self, request, session_slug):
        form = ConditionsForm(request.POST)
        print(request.GET)
        if form.is_valid():
            post = form.save(commit=False)
            post.session = Session.objects.get(identifier=session_slug)
            post.save()
            return redirect('analysis:session_detail', session_slug)
        return render(request, self.template_name, {'form':form})


class ConditionsUpdateView(UpdateView):
    template_name = 'analysis/conditions_form.html'
    form_class = ConditionsForm

    def get_object(self):
        conditions_pk = self.kwargs.get('conditions_pk')
        return get_object_or_404(Conditions, id=conditions_pk)

    def form_valid(self, form):
        session_slug = self.kwargs.get('session_slug')
        post = form.save(commit=False)
        post.session = Session.objects.get(identifier=session_slug)
        post.save()
        return redirect('analysis:session_detail', session_slug)


class ConditionsDeleteView(DeleteView):
    template_name = 'analysis/conditions_confirm_delete.html'
    context_object_name = 'condition'

    def get_object(self):
        # return HttpResponse('test')
        conditions_pk = self.kwargs.get('conditions_pk')
        return get_object_or_404(Conditions, pk=conditions_pk)

    def post(self, request, session_slug, conditions_pk):
        instance = get_object_or_404(Conditions, pk=conditions_pk)
        instance.delete()
        return redirect('analysis:session_detail', session_slug=session_slug)


# Samples
class SamplesListView(ListView):
    # context_object_name = 'conditions'
    model = models.Samples


class SamplesDetailView(DetailView):
    context_object_name = 'samples_detail'
    model = models.Samples
    template_name = 'analysis/samples_detail.html'


class SamplesCreateView(CreateView):
    template_name = 'analysis/samples_form.html'

    def get(self, request, session_slug, conditions_pk):
        form = SamplesForm
        context = {'form':form}
        session_pk = Session.objects.get(identifier=session_slug).id
        condition = Conditions.objects.get(pk=conditions_pk)
        con = condition.conditions
        print(f'\ncondition: {con}\n')
        print(f'replicates: {condition.no_replicates}')
        con_count = Samples.objects.filter(session_id=session_pk).filter(condition__conditions=con).count()
        print(f'\n{con_count}')
        if con_count >= condition.no_replicates:
            return redirect('analysis:session_detail', session_slug=session_slug)
        return render(request, self.template_name, context)

    def post(self, request, session_slug, conditions_pk):
        bound_form = SamplesForm(request.POST, request.FILES)
        # print(bound_form.is_valid())
        if bound_form.is_valid():
            # print(bound_form.cleaned_data)
            bound_post = bound_form.save(commit=False)
            bound_post.session = Session.objects.get(identifier=session_slug)
            bound_post.condition = Conditions.objects.get(pk=conditions_pk)
            bound_post.save()
            return redirect('analysis:session_detail', session_slug=session_slug)
        return render(request, self.template_name, {'form':form})


class SamplesUpdateView(UpdateView):
    template_name = 'analysis/samples_form.html'
    form_class = SamplesForm

    def get_object(self):
        samples_pk = self.kwargs.get('samples_pk')
        return get_object_or_404(Samples, pk=samples_pk)

    def form_valid(self, form):
        session_slug = self.kwargs.get('session_slug')
        form.save()
        return redirect('analysis:session_detail', session_slug=session_slug)


class SamplesDeleteView(DeleteView):
    template_name = 'analysis/samples_confirm_delete.html'
    context_object_name = 'sample'

    def get_object(self):
        samples_pk = self.kwargs.get('samples_pk')
        return get_object_or_404(Samples, pk=samples_pk)

    def post(self, request,  session_slug, samples_pk):
        instance = get_object_or_404(Samples, pk=samples_pk)
        instance.delete()
        return redirect('analysis:session_detail', session_slug=session_slug)


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

    def get(self, request, session_slug):
        form = WorkflowForm
        context = {'form':form}
        return render(request, self.template_name, context)

    def post(self, request, session_slug):
        form = WorkflowForm(request.POST)
        if form.is_valid():
            post = form.save(commit=False)
            post.session = Session.objects.get(identifier=session_slug)
            post.save()
            return redirect('analysis:session_detail', session_slug)
        return render(request, self.template_name, {'form':form})


class WorkflowUpdateView(UpdateView):
    template_name = 'analysis/workflow_form.html'
    form_class = WorkflowForm

    def get_object(self):
        workflow_pk = self.kwargs.get('workflow_pk')
        return get_object_or_404(Workflow, pk=workflow_pk)

    def form_valid(self, form):
        session_slug = self.kwargs.get('session_slug')
        form.save()
        return redirect('analysis:session_detail', session_slug=session_slug)


class WorkflowDeleteView(DeleteView):
    template_name = 'analysis/workflow_confirm_delete.html'
    context_object_name = 'workflow'

    def get_object(self):
        workflow_pk = self.kwargs.get('workflow_pk')
        return get_object_or_404(Workflow, pk=workflow_pk)

    def post(self, request,  session_slug, workflow_pk):
        instance = get_object_or_404(Workflow, pk=workflow_pk)
        instance.delete()
        return redirect('analysis:session_detail', session_slug=session_slug)
