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
from analysis.models import Session

class VisualizationIndexView(View):
    def get(self, request, session_slug):

        context = {'inject':'im the injection'}
        return render(request, 'visualization/index.html', context)

def WorkFlowOneView(request, session_slug, wf_key, *args, **kwargs):
     # 2cf7e0a9-c209-4dd8-a5ba-4b6bd07f6b91
    print(wf_key)
    # wf_dir = ['star_samtools_stringtie_prepde_deseq2', 'hisat2_samtools_stringtie_prepde_deseq2']
    # if wf_key == 1:
    #     wf_infile = wf_dir[0]
    # else:
    #     wf_infile = wf_dir[1]
    # session_output_dir = os.path.join(settings.DATA_DIR, session_slug, 'output', wf_infile)
    # print(session_output_dir)
    # session_output_csv = os.path.join(session_output_dir, os.listdir(session_output_dir)[0])
    # print(f'\n{session_output_csv}')

    session = Session.objects.get(identifier = session_slug)
    workflows = session.workflow_fk.all()
    wf_path = workflows[0].paths # will need to iterate through if multi wf selected
    wf_path = eval(wf_path) # EXTRACT DICT FROM STRING IF EXTRACTED FROM DB WITH QUOTES
    # print(f'\n{wf_path}')
    norm_path = wf_path['norm']
    DGE_path = wf_path['DGE'][0]
    print(f'\nNormal File Path {norm_path}')
    print(f'\nDGE File Path {DGE_path}')

    arr = []
    with open (session_output_csv) as csvFile:
        csvReader = csv.DictReader(csvFile)
        print(csvReader)
        for csvRow in csvReader:
            arr.append(csvRow)
        return JsonResponse(arr, safe=False)


# {"norm": "/project/data/rnaseq/Data/2cf7e0a9-c209-4dd8-a5ba-4b6bd07f6b91/output/salmonquant_salmoncount_deseq2/norm_count.csv", "DGE": ["/project/data/rnaseq/Data/2cf7e0a9-c209-4dd8-a5ba-4b6bd07f6b91/output//salmonquant_salmoncount_deseq2/Normal-Tumour_DGE_res.csv"]}

class DebugView(View):
    def get(self, request):
        return HttpResponse('success at the debug')


# https://stackoverflow.com/questions/26453916/passing-data-from-django-to-d3
# https://stackoverflow.com/questions/43617277/sending-json-data-from-django-to-d3-js
