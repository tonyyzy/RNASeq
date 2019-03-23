"""webportal URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.contrib import admin
from django.urls import include, path
from django.conf import settings
from django.conf.urls.static import static


urlpatterns = [
    path('admin/', admin.site.urls),
    # path('chaining', include('smart_selects.urls')),
    path('', include('analysis.urls')),
    # path('', include('pages.urls')),
    path('users/', include('users.urls')),
    path('visualization/<slug:session_slug>/', include('visualization.urls'))
]

if settings.DEBUG:
    urlpatterns += static(settings.IMG_URL, document_root=settings.IMG_ROOT)

# https://github.com/digi604/django-smart-selects
# https://stackoverflow.com/questions/29460078/how-to-use-django-smart-select
