from django.urls import path
from . import views


# urlpatterns = [
#     path('', views.index, name='index'),
#     # path('fastq/', views.get_fastq, name='fastq')
# ]

urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.home_view, name='home_view'),
    path('test/', views.test_view, name='test'),
    path('about/', views.about_view, name='about'),
    path('analysis/', views.analysis_view, name='analysis'),
    path('detail/', views.detail_view, name='detail'),
    path('db_view/', views.db_view, name='db'),
]
