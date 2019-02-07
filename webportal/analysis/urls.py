from django.urls import path, re_path
from . import views
# from analysis.views import analysisView
from django.views.generic import TemplateView

app_name = 'analysis'

urlpatterns = [
    # path('', TemplateView.as_view(template_name="about.html")),
    path('', views.home_view, name='home'),
    path('upload_session', views.session_view, name='upload_session'),
    # path('session_list', views.session_list_view, name='session_list'),
    # path('session_detail/<int:session_id>', views.session_detail_view, name='session_detail'),
    path('conditions_list', views.conditions_list_view, name='conditions_list'),
    # path('upload_conditions', views.conditions_view, name='upload_conditions'),
    path('upload_samples', views.samples_view, name='upload_samples'),
    path('upload_workflow', views.workflow_view, name='upload_workflow'),
    path('samples_list', views.samples_list_view, name='samples_list'),
    # path('conditions_list/<int:first>/', views.conditions_list_view, name='conditions_list'),
    path('conditions_list/<int:session_id>', views.conditions_list_view, name='conditions_list'),
    # re_path('^condition_list\/(?P<year>[0-9]{4})$', views.conditions_list_view,  name='conditions_list'),

    # url(r'^cbv$', views.cbv_view.as_view()),
    path('cbv', views.cbv_view.as_view()),
    path('session_list', views.SessionListView.as_view(), name='session_list'),
    path('session_detail/<int:pk>', views.SessionDetailView.as_view(), name='session_detail'),
    path('session_create', views.SessionCreateView.as_view(), name='session_create'),
    path('session_update/<int:pk>', views.SessionUpdateView.as_view(), name='session_update'),
    path('session_delete/<int:pk>', views.SessionDeleteView.as_view(), name='session_delete'),
]
# url(r'^update/(?P<pk>\d+)/$',views.SessionUpdateView.as_view(),name='update'),
