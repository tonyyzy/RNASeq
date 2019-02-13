from django.urls import path, re_path
from . import views
# from analysis.views import analysisView
from django.views.generic import TemplateView

app_name = 'analysis'

urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.SessionIndexView.as_view(), name='session_index'),
    path('session_detail/<int:pk>/conditions_update/<int:id>/', views.ConditionsUpdateView.as_view(), name='conditions_update'),


    path('session_list', views.SessionListView.as_view(), name='session_list'),
    path('session_detail/<int:pk>/', views.SessionDetailView.as_view(), name='session_detail'),
    path('session_create', views.SessionCreateView.as_view(), name='session_create'),
    path('session_update/<int:pk>/', views.SessionUpdateView.as_view(), name='session_update'),
    path('session_delete/<int:pk>/', views.SessionDeleteView.as_view(), name='session_delete'),

    path('conditions_list', views.ConditionsListView.as_view(), name='conditions_list'),
    path('conditions_detail/<int:pk>/', views.ConditionsDetailView.as_view(), name='conditions_detail'),
    path('conditions_create/<int:pk>/', views.ConditionsCreateView.as_view(), name='conditions_create'),
    # path('<int:pk>/conditions_update/<int:pk>/', views.ConditionsUpdateView.as_view(), name='conditions_update'),
    path('conditions_delete/<int:pk>/', views.ConditionsDeleteView.as_view(), name='conditions_delete'),

    path('samples_list', views.SamplesListView.as_view(), name='samples_list'),
    path('samples_detail/<int:pk>/', views.SamplesDetailView.as_view(), name='samples_detail'),
    path('samples_create/<int:pk>/', views.SamplesCreateView.as_view(), name='samples_create'),
    path('samples_update/<int:pk>/', views.SamplesUpdateView.as_view(), name='samples_update'),
    path('samples_delete/<int:pk>/', views.SamplesDeleteView.as_view(), name='samples_delete'),

    path('workflow_list', views.WorkflowListView.as_view(), name='workflow_list'),
    path('workflow_detail/<int:pk>/', views.WorkflowDetailView.as_view(), name='workflow_detail'),
    path('workflow_create/<int:pk>/', views.WorkflowCreateView.as_view(), name='workflow_create'),
    path('workflow_update/<int:pk>/', views.WorkflowUpdateView.as_view(), name='workflow_update'),
    path('workflow_delete/<int:pk>/', views.WorkflowDeleteView.as_view(), name='workflow_delete'),
]
