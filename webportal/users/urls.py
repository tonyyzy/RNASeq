from django.urls import path
from . import views



urlpatterns = [
    path('', views.regiser_view, name='regiser_view'),
    # path('test/', views.test_view, name='test'),
    ]
