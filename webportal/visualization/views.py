from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)
from django.http import JsonResponse
import csv , json
from django.conf import settings
import os

class VisualizationIndexView(View):
    def get(self, request, session_slug):
        context = {'inject':'im the injection'}
        return render(request, 'visualization/index.html', context)


def WorkFlowOneView(request, session_slug, wf_key, *args, **kwargs):
    # fcd72896-218a-45f8-be02-361b3c94e192
    print(wf_key)
    wf_dir = ['star_samtools_stringtie_prepde_deseq2', 'hisat2_samtools_stringtie_prepde_deseq2']
    if wf_key == 1:
        wf_infile = wf_dir[0]
    else:
        wf_infile = wf_dir[1]

    session_output_dir = os.path.join(settings.DATA_DIR, session_slug, 'output', wf_infile)
    print(session_output_dir)
    session_output_csv = os.path.join(session_output_dir, os.listdir(session_output_dir)[0])
    print(f'\n{session_output_csv}')
    arr = []
    with open (session_output_csv) as csvFile:
        csvReader = csv.DictReader(csvFile)
        print(csvReader)
        for csvRow in csvReader:
            arr.append(csvRow)
        return JsonResponse(arr, safe=False)



class DebugView(View):
    def get(self, request):
        return HttpResponse('success at the debug')


# https://stackoverflow.com/questions/26453916/passing-data-from-django-to-d3
# https://stackoverflow.com/questions/43617277/sending-json-data-from-django-to-d3-js
