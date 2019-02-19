from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from .forms import SessionSearchForm, SessionForm, WorkflowForm, SamplesForm, ConditionsForm
from analysis.models import Session, Workflow, Samples, Conditions
from . import models
from django.db.models import Q
# from django.views.generic import TemplateView
# from django.core.files.storage import FileSystemStorage
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)


# Session
class SessionIndexView(View):
    def get(self, request):
        # return HttpResponse('Session Index View')
        sessions = Session.objects.all()
        query = self.request.GET.get('q')
        if query:
            query = query.replace('-','')
            print(f'cleaned query: {query}')
            query_set_list = sessions.filter(
                            Q(identifier__icontains=query))

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
    # template_name = 'analysis/session_detail.html'
    model = models.Session

    def get_context_data(self, **kwargs):
        context = super(SessionDetailView, self).get_context_data(**kwargs)
        session_pk = self.kwargs.get('pk')
        samples = Samples.objects.all()
        session_samples = []
        for i in range(len(samples)):
            if (samples[i].condition.session.pk == session_pk):
                session_samples.append(samples[i])
        context['samples'] = session_samples
        # print(context)
        return context


class SessionCreateView(CreateView):
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

    def get(self, request, pk):
        form = ConditionsForm
        context = {'form':form}
        return render(request, self.template_name, context)

    def post(self, request, pk):
        form = ConditionsForm(request.POST)
        print(request.GET)
        if form.is_valid():
            post = form.save(commit=False)
            sessions = Session.objects.all()
            adjusted_pk = self.kwargs.get('pk')-1 # database object starts index from 1 not 0
            # print(f'\n{sessions[adjusted_pk]}') # returns a session object NOT a string
            post.session = sessions[adjusted_pk]
            post.save()
            return redirect('analysis:session_detail', pk)
        return render(request, self.template_name, {'form':form})


class ConditionsUpdateView(UpdateView):
    template_name = 'analysis/conditions_form.html'
    form_class = ConditionsForm

    def get_object(self):
        conditions_pk = self.kwargs.get('conditions_pk')
        return get_object_or_404(Conditions, id=conditions_pk)

    def form_valid(self, form):
        session_pk = self.kwargs.get('session_pk')
        post = form.save(commit=False)
        post.session = Session.objects.get(pk=session_pk)
        post.save()
        return redirect('analysis:session_detail', pk=session_pk)


class ConditionsDeleteView(DeleteView):
    template_name = 'analysis/conditions_confirm_delete.html'
    context_object_name = 'condition'

    def get_object(self):
        conditions_pk = self.kwargs.get('conditions_pk')
        return get_object_or_404(Conditions, id=conditions_pk)

    def post(self, request, session_pk, conditions_pk):
        instance = get_object_or_404(Conditions, id=conditions_pk)
        instance.delete()
        return redirect('analysis:session_detail', pk=session_pk)


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

    def get(self, request, session_pk, conditions_pk):
        form = SamplesForm
        context = {'form':form}
        return render(request, self.template_name, context)

    def post(self, request, session_pk, conditions_pk):
        form = SamplesForm
        bound_form = SamplesForm(request.POST, request.FILES)
        # print(bound_form.is_valid())
        if bound_form.is_valid():
            print(bound_form.cleaned_data)
            session = Session.objects.get(pk=session_pk)
            condition = Conditions.objects.get(pk=conditions_pk)
            bound_post = bound_form.save(commit=False)
            bound_post.session = session
            bound_post.condition = condition
            bound_post.save()
            return redirect('analysis:session_detail', pk=session_pk)
        return render(request, self.template_name, {'form':form})


class SamplesUpdateView(UpdateView):
    template_name = 'analysis/samples_form.html'
    form_class = SamplesForm

    def get_object(self):
        samples_pk = self.kwargs.get('samples_pk')
        return get_object_or_404(Samples, pk=samples_pk)

    def form_valid(self, form):
        session_pk = self.kwargs.get('session_pk')
        post = form.save(commit=False)
        # post.session = Session.objects.get(pk=session_pk)
        post.save()
        return redirect('analysis:session_detail', pk=session_pk)


class SamplesDeleteView(DeleteView):
    template_name = 'analysis/samples_confirm_delete.html'
    context_object_name = 'sample'

    def get_object(self):
        samples_pk = self.kwargs.get('samples_pk')
        return get_object_or_404(Samples, pk=samples_pk)

    def post(self, request,  session_pk, samples_pk):
        instance = get_object_or_404(Samples, pk=samples_pk)
        instance.delete()
        return redirect('analysis:session_detail', pk=session_pk)


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

    def get(self, request, session_pk):
        form = WorkflowForm
        context = {'form':form}
        return render(request, self.template_name, context)

    def post(self, request, session_pk):
        form = WorkflowForm(request.POST)
        if form.is_valid():
            session = Session.objects.filter(id=session_pk)
            post = form.save(commit=False)
            post.session = Session.objects.filter(id=session_pk)[0]
            post.save()
            return redirect('analysis:session_detail', session_pk)
        return render(request, self.template_name, {'form':form})


class WorkflowUpdateView(UpdateView):
    template_name = 'analysis/workflow_form.html'
    form_class = WorkflowForm

    def get_object(self):
        # return HttpResponse('hi there old sport')
        workflow_pk = self.kwargs.get('workflow_pk')
        return get_object_or_404(Workflow, pk=workflow_pk)

    def form_valid(self, form):
        session_pk = self.kwargs.get('session_pk')
        post = form.save(commit=False)
        # post.session = Session.objects.get(pk=session_pk)
        post.save()
        return redirect('analysis:session_detail', pk=session_pk)


class WorkflowDeleteView(DeleteView):
    template_name = 'analysis/workflow_confirm_delete.html'
    context_object_name = 'workflow'

    def get_object(self):
        workflow_pk = self.kwargs.get('workflow_pk')
        return get_object_or_404(Workflow, pk=workflow_pk)

    def post(self, request,  session_pk, workflow_pk):
        instance = get_object_or_404(Workflow, pk=workflow_pk)
        instance.delete()
        return redirect('analysis:session_detail', pk=session_pk)
