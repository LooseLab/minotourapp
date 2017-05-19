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
from reads.models import FastqRead
from reads.models import RunStatistic
from reads.models import FastqReadType


class FastqReadSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('id', 'run_id', 'read_id', 'read', 'channel', 'barcode', 'sequence', 'quality', 'is_pass',
                  'start_time', 'type')
        read_only = ('id',)


class RunStatisticSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = RunStatistic
        fields = ('run_id', 'sample_time', 'total_length', 'max_length', 'min_length', 'average_length',
                  'number_of_reads', 'number_of_channels', 'type')


class FastqReadTypeSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqReadType
        fields = ('id', 'name',)
        read_only = ('id',)


class RunStatisticViewSet(viewsets.ModelViewSet):
    queryset = RunStatistic.objects.all()
    serializer_class = RunStatisticSerializer


class FastqReadViewSet(viewsets.ModelViewSet):
    queryset = FastqRead.objects.all()
    serializer_class = FastqReadSerializer


class FastqReadTypeViewSet(viewsets.ModelViewSet):
    queryset = FastqReadType.objects.all()
    serializer_class = FastqReadTypeSerializer


router = routers.DefaultRouter()
router.register(r'reads', FastqReadViewSet)
router.register(r'readtypes', FastqReadTypeViewSet)
router.register(r'statistics', RunStatisticViewSet)


urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^api/v1/', include(router.urls)),
    url(r'^', include('reads.urls')),
]
