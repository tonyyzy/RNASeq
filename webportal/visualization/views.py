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
        session = Session.objects.get(identifier = session_slug)
        workflows = session.workflow_fk.all()

        # workflow_detail = [workflows[0]]
        context = {'session': session, 'workflows': workflows}
        return render(request, 'visualization/index.html', context)

def WorkFlowOneView(request, session_slug, workflow_slug):
    session = Session.objects.get(identifier = session_slug)
    workflows = session.workflow_fk.all()
    selected_wf_list = workflow_slug[:-1].split("_")
    selected_wf = {workflows[int(i)-1] for i in selected_wf_list} # adds user selected workflow OBJECTS(s) to list
    selected_DGE = []
    for wf in selected_wf:
        wf_path = wf.paths
        wf_path = eval(wf_path) # only required if workflow path has been manually added to database in which case it is stored as a string
        wf_csv = wf_path['DGE'][0] # pulling out first here. Would need to loop if more than two sample conditions are defined.
        selected_DGE.append(wf_csv)

    arr = []
    for index, DGE_csv in enumerate(selected_DGE):
        print(f'\n{DGE_csv}')
        with open (DGE_csv) as in_file:
            csvReader = csv.DictReader(in_file)
            # print(csvReader)
            for csvRow in csvReader:
                csvRow['dataset_index'] = index
                print(csvRow)
                # print(csvRow)
                arr.append(csvRow)

    # print(arr)
    return JsonResponse(arr, safe=False)

# all_workflows = []
# for i in workflows:
#     if i.paths == '':
#         all_workflows.append('')
#     all_workflows.append(json.loads(i.paths))
# all_workflows = [all_workflows[i] for i in range(length(all_workflows)) if i in selected_wf]
# print(all_workflows)


# 128bf2d4-a9d1-4cf4-91ed-c897099467f8
# {"norm": "/project/home18/ph2417/gitWorkspace/Data/128bf2d4-a9d1-4cf4-91ed-c897099467f8/hisat_samtools_stringtie_prepde_deseq/norm_count.csv", "DGE": ["/project/home18/ph2417/gitWorkspace/Data/128bf2d4-a9d1-4cf4-91ed-c897099467f8/hisat_samtools_stringtie_prepde_deseq/DGE_results_1.csv"]}


class DebugView(View):
    def get(self, request, session_slug, workflow_slug):
        return HttpResponse(workflow_slug)


# https://stackoverflow.com/questions/26453916/passing-data-from-django-to-d3
# https://stackoverflow.com/questions/43617277/sending-json-data-from-django-to-d3-js
