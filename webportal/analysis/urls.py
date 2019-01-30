from django.urls import path
from . import views
# from analysis.views import analysisView
from django.views.generic import TemplateView

# urlpatterns = [
#     path('', views.index, name='index'),
#     # path('fastq/', views.get_fastq, name='fastq')
# ]

urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.index, name='index'),
    # path('test', views.test, name='test'),
    path('file_submit/', views.get_files, name='file_submit'),
    path('view_db/', views.my_view, name='my_view'),
]
