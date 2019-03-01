from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from .forms import SessionSearchForm, SessionForm, SessionSubmitForm, WorkflowForm, SamplesForm, ConditionsForm, DebugForm, GenomeForm
from analysis.models import Session, Workflow, Samples, Conditions, The_Debug
from . import models
from django.db.models import Q
from django.contrib import messages
import os
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
                print(f'\n{session_query}')
                messages.success(request, 'success')
                context = {'session_query':session_query}
                return render(request, 'analysis/index.html', context)
            except:
                messages.warning(request, 'Session ID not found')
                context = {'invalid_query': query}
                return render(request, 'analysis/index.html', context)
        return render(request, 'analysis/index.html')


class GenomeCreateView(CreateView):
    template_name = 'analysis/genome_form.html'

    def get(self, request):
        form = GenomeForm
        # return HttpResponse('success')
        return render(request, self.template_name, {'form':form})

    def post(self, request):
        form = GenomeForm
        bound_form = GenomeForm(request.POST, request.FILES)
        print(bound_form)
        # print(bound_form.cleaned_data)
        print(bound_form.is_valid())
        if bound_form.is_valid():
            post = bound_form.save()
            print(f'\n{post}')
            post.save()
            return redirect('analysis:session_index')
        return render(request, self.template_name, {'form':form})


class SessionListView(ListView):
    # context_object_name = 'session'
    # by default context object used to access database object within template is is lower(model_name)_list
    model = models.Session


class SessionDetailView(View):
    template_name = 'analysis/session_detail.html'

    def get(self, request, session_slug): # loads the session_detail template with the selected session object loaded as 'instance' and upload associated with that instance loaded from 'form'
        instance = get_object_or_404(Session, identifier=session_slug)
        form = SessionSubmitForm(request.POST or None, instance = instance)
        try:
            session = Session.objects.get(identifier=session_slug)
        except session.DoesNotExist:
            raise Http404('Session does not exist')

        context = {'session_detail':session, 'key2':'val2', 'form':form}
        return render(request, self.template_name, context)


    def post(self, request, session_slug):
        instance = get_object_or_404(Session, identifier=session_slug)
        bound_form = SessionSubmitForm(request.POST or None, instance = instance)
        if bound_form.is_valid():
            post = bound_form.save(commit=False)
            post.status = True
            post.save()
            messages.success(request, 'Upload Successful')
            print(f'\n{post}')
            return redirect('analysis:session_detail', session_slug)
        return redirect('analysis:session_detail', session_slug)


class SessionCreateView(CreateView):
    template_name = 'analysis/session_form.html'

    def get(self, request):
        form = SessionForm
        return render(request, self.template_name, {'form':form})

    def post(self, request):
        form = SessionForm
        bound_form = SessionForm(request.POST, request.FILES)
        if bound_form.is_valid():
            post = bound_form.save()
            return redirect('analysis:session_detail', session_slug=post.identifier)
        return render(request, self.template_name, {'form':form})


class SessionUpdateView(UpdateView):
    # fields = ('organism','genome','fasta_file', 'annotation_file',)
    # model = models.Session

    template_name = 'analysis/session_form.html'
    form_class = SessionForm

    def get_object(self):
        session_slug = self.kwargs.get('session_slug')
        print('object retreived success object retreived success  object retreived success ')
        return get_object_or_404(Session, identifier=session_slug)

    def form_valid(self, form):
        session_slug = self.kwargs.get('session_slug')
        post = form.save(commit=False)
        post.save()
        return redirect('analysis:session_detail', session_slug)



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
        condition_obj = Conditions.objects.get(pk=conditions_pk)
        row_condition = condition_obj.conditions
        print(f'\ncondition: {row_condition}\n')
        print(f'replicates: {condition_obj.no_replicates}')
        row_condition_count = Samples.objects.filter(session_id=session_pk).filter(condition__conditions=row_condition).count()
        if row_condition_count >= condition_obj.no_replicates:
            messages.warning(request, f'Error: {condition_obj.no_replicates} {condition_obj.conditions} sample(s) already uploaded.')
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


class DebugView(CreateView):

    def get(self, request):
        form = DebugForm
        debugs = The_Debug.objects.all()
        current_path = os.getcwd()
        onlyfiles = [f for f in os.listdir(current_path) if os.path.isfile(os.path.join(current_path, f))]
        return render(request, 'analysis/debug_page.html', {'form':form, 'debugs':debugs, 'onlyfiles':onlyfiles})

    def post(self, request):
        bound_form = DebugForm(request.POST or None, request.FILES)
        print(f'\n{bound_form.is_valid()}')
        if bound_form.is_valid():
            bound_form.save()
            print(f'\n{bound_form}')
            messages.success(request, 'Upload Successful')
            # print(f'\n{post}')
            return redirect('analysis:debug_view')
        return redirect('analysis:debug_view')


#
# class SessionCreateView(CreateView):
#     template_name = 'analysis/session_form.html'
#
#     def get(self, request):
#         form = SessionForm
#         return render(request, self.template_name, {'form':form})
#
#     def post(self, request):
#         form = SessionForm
#         bound_form = SessionForm(request.POST, request.FILES)
#         if bound_form.is_valid():
#             post = bound_form.save()
#             return redirect('analysis:session_detail', session_slug=post.identifier)
#         return render(request, self.template_name, {'form':form})
