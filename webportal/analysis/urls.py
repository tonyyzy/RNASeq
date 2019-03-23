from django.urls import path, re_path
from . import views
from django.views.generic import TemplateView

app_name = 'analysis'

urlpatterns = [

    path('', views.SessionIndexView.as_view(), name='session_index'),
    path('debug/', views.DebugView.as_view(), name='debug_view'),
    path('genome_create/', views.GenomeCreateView.as_view(), name='genome_view'),

    # path('session_detail/<slug:session_slug>/workflow_update/<int:workflow_pk>/filtered_analysis/<slug:assembler_slug>/', views.filterAnalysis, name='filtered_analysis'),

    path('session_list', views.SessionListView.as_view(), name='session_list'),
    path('session_detail/<slug:session_slug>/', views.SessionDetailView.as_view(), name='session_detail'),
    path('session_detail/<slug:session_slug>/svg_download', views.SVGDownload, name='svg_download'),

    path('session_create', views.SessionCreateView.as_view(), name='session_create'),
    path('session_update/<slug:session_slug>/', views.SessionUpdateView.as_view(), name='session_update'),
    path('session_delete/<int:pk>/', views.SessionDeleteView.as_view(), name='session_delete'),

    path('condition_list', views.ConditionListView.as_view(), name='condition_list'),
    path('conditions_detail/<int:pk>/', views.ConditionsDetailView.as_view(), name='conditions_detail'),
    path('session_detail/<slug:session_slug>/conditions_create/', views.ConditionsCreateView.as_view(), name='conditions_create'),
    path('session_detail/<slug:session_slug>/conditions_update/<int:conditions_pk>/', views.ConditionsUpdateView.as_view(), name='conditions_update'),
    path('session_detail/<slug:session_slug>/conditions_delete/<int:conditions_pk>/', views.ConditionsDeleteView.as_view(), name='conditions_delete'),

    path('samples_list', views.SamplesListView.as_view(), name='samples_list'),
    path('samples_detail/<int:pk>/', views.SamplesDetailView.as_view(), name='samples_detail'),
    path('session_detail/<slug:session_slug>/samples_create/<int:conditions_pk>/', views.SamplesCreateView.as_view(), name='samples_create'),
    path('session_detail/<slug:session_slug>/samples_update/<int:samples_pk>/', views.SamplesUpdateView.as_view(), name='samples_update'),
    path('session_detail/<slug:session_slug>/samples_delete/<int:samples_pk>/', views.SamplesDeleteView.as_view(), name='samples_delete'),

    path('workflow_list', views.WorkflowListView.as_view(), name='workflow_list'),
    path('workflow_detail/<int:pk>/', views.WorkflowDetailView.as_view(), name='workflow_detail'),
    path('session_detail/<slug:session_slug>/workflow_create/', views.WorkflowCreateView.as_view(), name='workflow_create'),
    path('session_detail/<slug:session_slug>/workflow_create/filtered_assembler/<slug:mapper_slug>', views.filterAssembler, name='filtered_assembler'),
    path('session_detail/<slug:session_slug>/workflow_create/filtered_analysis/<slug:assembler_slug>', views.filterAnalysis, name='filtered_analysis'),

    path('session_detail/<slug:session_slug>/workflow_update/<int:workflow_pk>/', views.WorkflowUpdateView.as_view(), name='workflow_update'),
    path('session_detail/<slug:session_slug>/workflow_update/<int:workflow_pk>/filtered_assembler/<slug:mapper_slug>', views.filterAssemblerUpdate, name='filtered_assembler'),
    path('session_detail/<slug:session_slug>/workflow_update/<int:workflow_pk>/filtered_analysis/<slug:assembler_slug>', views.filterAnalysisUpdate, name='filtered_analysis'),

    path('session_detail/<slug:session_slug>/workflow_delete/<int:workflow_pk>/', views.WorkflowDeleteView.as_view(), name='workflow_delete'),

]
# 04a2d270-2d00-4370-8538-08dc33fbba6e
