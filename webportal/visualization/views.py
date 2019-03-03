from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)


class VisualizationIndexView(View):
    def get(self, request):
        # return HttpResponse('success at the viz')
        return render(request, 'visualization/index.html', {'inject':'im the injection'})


class DebugView(View):
    def get(self, request):
        return HttpResponse('success at the debug')
