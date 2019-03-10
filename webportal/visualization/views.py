from django.shortcuts import render, redirect, get_object_or_404
from django.http import HttpResponse
from django.urls import reverse, reverse_lazy
from django.views.generic import (View,TemplateView,
                                ListView,DetailView,
                                CreateView,DeleteView,
                                UpdateView)
from django.http import JsonResponse
import csv , json


class VisualizationIndexView(View):
    def get(self, request):
        context = {'inject':'im the injection'}
        return render(request, 'visualization/index.html', context)


def GetDataView(request, *args, **kwargs):
    json_data = {"greeting": "world", "foo": "bar", "gatsby": "hello old sport"}
    csvFilePath = "/home/patrick/Code/GitWorkSpace/myApp/testPlot/DGE_results_1_head.csv" # will convert to use session specifi path once d3.json(data) is working
    # jsonFilePath = "j_test.json"
    arr = []
    #read the csv and add the arr to a array
    with open (csvFilePath) as csvFile:
        csvReader = csv.DictReader(csvFile)
        print(csvReader)
        for csvRow in csvReader:
            arr.append(csvRow)
    print(arr)
    # return JsonResponse(json_data)
    return JsonResponse(arr, safe=False)


class DebugView(View):
    def get(self, request):
        return HttpResponse('success at the debug')


# https://stackoverflow.com/questions/26453916/passing-data-from-django-to-d3
# https://stackoverflow.com/questions/43617277/sending-json-data-from-django-to-d3-js
