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
import pandas
from collections import OrderedDict, defaultdict

class VisualizationIndexView(View):
    def get(self, request, session_slug):
        session = Session.objects.get(identifier = session_slug)
        workflows = session.workflow_fk.all()
        # workflow = workflows[0]

        wf_dict = {}
        labels = []
        for wf in workflows:
            labels.append(wf.label)
            print(f'\n workflow analysis: {wf.analysis}')
            if wf.analysis == 'dexseq':
                path = eval(wf.paths)
                DGE_list = path['DEE']

                DGE_csv_list = []
                for DGE_csv in DGE_list:
                    csv_name = DGE_csv.split('/')[-1]
                    DGE_csv_list.append(csv_name)
                wf_dict[wf.label] = DGE_csv_list
            else:
                path = eval(wf.paths)
                print(path)
                DEE_list = path['DGE']
                DEE_csv_list = []
                for DEE_csv in DEE_list:
                    csv_name = DEE_csv.split('/')[-1]
                    DEE_csv_list.append(csv_name)
                wf_dict[wf.label] = DEE_csv_list


        # print(wf_dict)
        # for workflow in workflows:
            # path = workflow.paths
        # workflow_detail = [workflows[0]]
        context = {'session': session, 'workflows': workflows, 'labels':labels, 'wf_dict': wf_dict}
        return render(request, 'visualization/index.html', context)

def WorkflowData(request, session_slug, workflow_slug):
    session = Session.objects.get(identifier = session_slug)
    workflows = session.workflow_fk.all()
    selected_wf_list = workflow_slug[:-1].split("_")
    selected_wf = {workflows[int(i)-1] for i in selected_wf_list} # adds user selected workflow OBJECTS(s) to list
    print(f'\n{selected_wf}')
    selected_DGE = []

    for wf in selected_wf:
        wf_path = wf.paths
        wf_path = eval(wf_path)
        print(f'\n{wf_path}')

        print(f'workflow analysis: {wf.analysis}')
        if wf.analysis == 'deseq2':
            DGE_csv_list = wf_path['DGE']

            print(f'\DGE_csv_list: {DGE_csv_list}')
            for DGE_csv in DGE_csv_list:
                selected_DGE.append(DGE_csv)

        if wf.analysis == 'dexseq':
            pass

    arr = []
    for index, DGE_csv in enumerate(selected_DGE):
        print(f'\n{DGE_csv}')
        with open (DGE_csv) as in_file:
            csvReader = pandas.read_csv(in_file)
            csvReader['dataset_index'] = index
            arr.append(csvReader)
            # print(arr)
            # csvReader = csv.DictReader(in_file)
            # print(csvReader)
            # for csvRow in csvReader:
            #     csvRow['dataset_index'] = index
            #     print(csvRow)
            #     arr.append(csvRow)
    # print(arr)
    final = pandas.concat(arr, keys=range(len(arr)))
    # print(final)
    # final = final.round(8)
    final = final.to_json(orient='records')
    return HttpResponse(final)
    # return render(final, safe=False)
    # print(arr)
    # return JsonResponse(final, safe=False)



# 128bf2d4-a9d1-4cf4-91ed-c897099467f8
# {"norm": "/project/home18/ph2417/gitWorkspace/Data/128bf2d4-a9d1-4cf4-91ed-c897099467f8/hisat_samtools_cufflinks_prepde_deseq/DGE_results_1.csv", "DGE": ["/project/home18/ph2417/gitWorkspace/Data/128bf2d4-a9d1-4cf4-91ed-c897099467f8/hisat_samtools_cufflinks_prepde_deseq/DGE_results_1.csv"]}

# {"norm": "/project/home18/ph2417/gitWorkspace/Data/128bf2d4-a9d1-4cf4-91ed-c897099467f8/star_samtools_stringtie_prepde_deseq/norm_count.csv", "DGE": ["/project/home18/ph2417/gitWorkspace/Data/128bf2d4-a9d1-4cf4-91ed-c897099467f8/star_samtools_stringtie_prepde_deseq/DGE_results_1.csv"]}






class DebugView(View):
    # def get(self, request, session_slug, workflow_slug):
    #     return HttpResponse(workflow_slug)

    def get(request, session_slug, workflow_slug):
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


# https://stackoverflow.com/questions/26453916/passing-data-from-django-to-d3
# https://stackoverflow.com/questions/43617277/sending-json-data-from-django-to-d3-js
