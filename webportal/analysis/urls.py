from django.urls import path
from . import views
# from analysis.views import analysisView
from django.views.generic import TemplateView


urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.test_view, name='test_view'),
    path('upload', views.samples_view, name='upload'),
    path('file_submit/', views.get_files, name='file_submit'),
    path('view_db/', views.my_view, name='my_view'),
]
