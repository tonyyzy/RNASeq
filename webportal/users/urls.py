from django.urls import path
from . import views


app_name = 'users'

urlpatterns = [
    path('signup', views.regiser_view, name='regiser_view'),
    path('login', views.login_view, name='login_view'),
    path('logout', views.logout_view, name='logout_view'),
    # path('', views.test_view, name='test'),
    ]
