from django.urls import path, re_path
from . import views
from django.views.generic import TemplateView

app_name = 'visualization'

urlpatterns = [
    path('', views.VisualizationIndexView.as_view(), name='visualization_index'),
    # path('wf_data/<slug:workflow_slug>/', views.WorkflowData, name='wf_data'),cd
    path('wf_data_mod/<slug:workflow_slug>/', views.WorkflowDataMod, name='wf_data_mod'),
    path('wf_file_download/<slug:workflow_slug>/', views.wfDownload, name='wf_file_download'),
    path('debug/', views.debugFunction, name='debug_func'),
    # path('wf_one/<slug:workflow_slug>/', views.WorkFlowOneView, name='wf_one'),
    # path('wf_two/', views.WorkFlowTwoView, name='wf_two'),
]
