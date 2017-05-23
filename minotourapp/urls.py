"""minotourapp URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import include
from django.conf.urls import url
from django.contrib import admin
from rest_framework import routers
from rest_framework import serializers
from rest_framework import viewsets

from reads.models import RunStatistic


class RunStatisticSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = RunStatistic
        fields = ('run_id', 'sample_time', 'total_length', 'max_length', 'min_length',
                  'read_count', 'type')


class RunStatisticViewSet(viewsets.ModelViewSet):
    queryset = RunStatistic.objects.all()
    serializer_class = RunStatisticSerializer


router = routers.DefaultRouter()
router.register(r'statistics', RunStatisticViewSet)

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^api/v1/', include(router.urls)),
    url(r'^web/', include('web.urls')),
    url(r'^', include('reads.urls')),
]
