from django.urls import path
from . import views
# from analysis.views import analysisView
from django.views.generic import TemplateView


urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.test_view, name='test_view'),
    path('upload_session', views.session_view, name='upload'),
    path('upload_samples', views.samples_view, name='upload'),
    path('upload_workflow', views.workflow_view, name='upload'),
]
