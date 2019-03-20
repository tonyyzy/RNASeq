from django.urls import path, re_path
from . import views
from django.views.generic import TemplateView

app_name = 'visualization'

urlpatterns = [
    path('', views.VisualizationIndexView.as_view(), name='visualization_index'),
    path('debug/<slug:workflow_slug>/', views.DebugView.as_view(), name='debug_view'),
    path('wf_data/<slug:workflow_slug>/', views.WorkflowData, name='wf_data'),
    # path('wf_one/<slug:workflow_slug>/', views.WorkFlowOneView, name='wf_one'),
    # path('wf_two/', views.WorkFlowTwoView, name='wf_two'),
]
#test
