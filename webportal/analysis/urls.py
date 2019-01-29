from django.urls import path
from . import views

# urlpatterns = [
#     path('', views.index, name='index'),
#     # path('fastq/', views.get_fastq, name='fastq')
# ]

urlpatterns = [
    path('', views.index, name='index'),
    path('test', views.test, name='test'),
    path('file_submit/', views.get_files, name='file_submit'),
    # path('session/', views.session_form_upload, name='session')
]
