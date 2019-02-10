from django.urls import path, re_path
from . import views
# from analysis.views import analysisView
from django.views.generic import TemplateView

app_name = 'analysis'

urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.SessionListView.as_view(), name='session_list'),

    # function based views
    # path('upload_conditions', views.conditions_view, name='upload_conditions'),
    # path('upload_samples', views.samples_view, name='upload_samples'),
    # path('upload_workflow', views.workflow_view, name='upload_workflow'),
    # path('samples_list', views.samples_list_view, name='samples_list'),


    # class based views
    # url(r'^cbv$', views.cbv_view.as_view()),
    path('cbv', views.cbv_view.as_view()),
    path('session_list', views.SessionListView.as_view(), name='session_list'),
    path('session_detail/<int:pk>/', views.SessionDetailView.as_view(), name='session_detail'),
    path('session_create', views.SessionCreateView.as_view(), name='session_create'),
    path('session_update/<int:pk>/', views.SessionUpdateView.as_view(), name='session_update'),
    path('session_delete/<int:pk>/', views.SessionDeleteView.as_view(), name='session_delete'),

    path('conditions_list', views.ConditionsListView.as_view(), name='conditions_list'),
    path('conditions_detail/<int:pk>/', views.ConditionsDetailView.as_view(), name='conditions_detail'),
    path('conditions_create/<int:pk>/', views.ConditionsCreateView.as_view(), name='conditions_create'),
    path('conditions_update/<int:pk>/', views.ConditionsUpdateView.as_view(), name='conditions_update'),
    path('conditions_delete/<int:pk>/', views.ConditionsDeleteView.as_view(), name='conditions_delete'),

    path('samples_list', views.SamplesListView.as_view(), name='samples_list'),
    path('samples_detail/<int:pk>/', views.SamplesDetailView.as_view(), name='samples_detail'),
    path('samples_create/<int:pk>/', views.SamplesCreateView.as_view(), name='samples_create'),
    path('samples_update/<int:pk>/', views.SamplesUpdateView.as_view(), name='samples_update'),
    path('samples_delete/<int:pk>/', views.SamplesDeleteView.as_view(), name='samples_delete'),
]
