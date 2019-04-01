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
from analysis.models import Session, Workflow
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
                DEE_list = path['DEE']
                DEE_csv_list = []
                for DEE_csv in DEE_list:
                    csv_name = DEE_csv.split('/')[-1]
                    DEE_csv_list.append(csv_name)
                wf_dict[wf] = DEE_csv_list

            else:
                try:
                    path = eval(wf.paths)
                    DGE_list = path['DGE']
                    DGE_csv_list = []
                    for DGE_csv in DGE_list:
                        csv_name = DGE_csv.split('/')[-1]
                        DGE_csv_list.append(csv_name)
                    wf_dict[wf] = DGE_csv_list
                except:
                    pass

        # print(wf_dict)
        # for workflow in workflows:
            # path = workflow.paths
        # workflow_detail = [workflows[0]]
        context = {'session': session, 'workflows': workflows, 'labels':labels, 'wf_dict': wf_dict}
        return render(request, 'visualization/index.html', context)

def WorkflowDataMod(request, session_slug, workflow_slug):
    session = Session.objects.get(identifier = session_slug)
    workflows = Workflow.objects.all()
    print(f'\nworkflow slug: {workflow_slug}')
    # print(f'\n{workflows.values}')
    selected_wf_list = workflow_slug[:-1].split('_')
    # selected_wf_list = [wf for wf in selected_wf_list if wf % 2 == 1]

    selected_wf_list = selected_wf_list[1::2]
    print(f'\nselected workflow list: {selected_wf_list}')
    selected_wf = []
    selected_file = {}
    for wf in selected_wf_list:
        wf, wf_file = wf.split('-')

        selected_wf.append(wf)
        if wf in selected_file.keys():
            selected_file[wf].append(wf_file)
        else:
            selected_file[wf] = []
            selected_file[wf].append(wf_file)

    print(f'\nselected wf: {selected_wf}')
    print(f'\nselected file: {selected_file}')

    selected_csv = []
    selected_wf = list(set(selected_wf))

    for i in selected_wf:
        workflow = workflows.get(id=i)
        # print(f'\nTHE WORKFLOW: {workflow}')
        paths = workflow.paths
        paths = eval(paths)
        paths = paths['DGE']
        # print(paths)
        # print(paths[0])
        # print(i)
        # print([selected_file[i]])
        for e in selected_file[i]:
            print(f"this is test {int(e) - 1}")
            selected_csv.append(paths[int(e) - 1])

    # selected_csv = selected_csv[0]
    # print(f'\n the selected csv: {selected_csv}')

    arr = []
    for index, DGE_csv in enumerate(selected_csv):
        print('test')

        print(f'\nthe DGE_csv: {DGE_csv}')
        with open (DGE_csv) as in_file:
            csvReader = pandas.read_csv(in_file)
            csvReader['dataset_index'] = index
            arr.append(csvReader)
    final = pandas.concat(arr, keys=range(len(arr)))
    print(final.shape)
    final = final.to_json(orient='records')
    return HttpResponse(final)
    # return JsonResponse(arr, safe=False)


def wfDownload(request, session_slug, workflow_slug):
        session = Session.objects.get(identifier = session_slug)
        workflows = Workflow.objects.all()
        return HttpResponse(workflow_slug)
        print(f'\nworkflow slug: {workflow_slug}')
        # print(f'\n{workflows.values}')
        selected_wf_list = workflow_slug[:-1].split('_')
        # selected_wf_list = [wf for wf in selected_wf_list if wf % 2 == 1]
        selected_wf_list = selected_wf_list[1::2]
        print(f'\n{selected_wf_list}')
    # print(f'\n SVG Downlod called')
    # img_path = os.path.join(settings.BASE_DIR, 'static/images', session_slug, 'workflow.svg')
    # img_wrapper = FileWrapper(open(img_path,'rb'))
    # response = HttpResponse(img_wrapper)
    # response['X-Sendfile'] = img_path
    # response['Content-Length'] = os.stat(img_path).st_size
    # response['Content-Disposition'] = 'attachment; filename=workflow.svg'
    # return response

#     with open (DGE_csv) as in_file:
#         csvReader = pandas.read_csv(in_file)
#         csvReader['dataset_index'] = index
#         arr.append(csvReader)
# final = pandas.concat(arr, keys=range(len(arr)))
# final = final.to_json(orient='records')
# return HttpResponse(final)


def debugFunction(request, session_slug):
    selected_csv = ['/project/data/rnaseq/Data/4e4dd873-f4a8-45a6-8a74-7b67173568b4/output//star_samtools_stringtie_prepde_deseq2/SY5Y37-6C37_DGE_res.csv']
    print(f'\n the selected csv: {selected_csv}')
    arr = []
    for index, DGE_csv in enumerate(selected_csv):
        print(f'\nthe DGE_csv: {DGE_csv}')
        with open (DGE_csv) as in_file:
            csvReader = csv.DictReader(in_file)
            for csvRow in csvReader:
                arr.append(csvRow)
    return JsonResponse(arr, safe=False)





# https://stackoverflow.com/questions/26453916/passing-data-from-django-to-d3
# https://stackoverflow.com/questions/43617277/sending-json-data-from-django-to-d3-js
