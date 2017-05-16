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

from rest_framework import routers, serializers, viewsets
from reads.models import FastqRead, RunStatistic
from reads.models import MinionRun


class MinionRunSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = MinionRun
        fields = ('run_name', 'run_id', 'barcode')


class FastqReadSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = FastqRead
        fields = ('run_id', 'read_id', 'read', 'channel', 'barcode', 'sequence', 'quality', 'status', 'start_time')


class RunStatisticSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = RunStatistic
        fields = ('run_id', 'sample_time', 'total_length', 'max_length', 'min_length', 'average_length', 'number_of_reads', 'number_of_channels') #, 'type')


class RunStatisticViewSet(viewsets.ModelViewSet):
    queryset = RunStatistic.objects.all()
    serializer_class = RunStatisticSerializer


class MinionRunViewSet(viewsets.ModelViewSet):
    queryset = MinionRun.objects.all()
    serializer_class = MinionRunSerializer


class FastqReadViewSet(viewsets.ModelViewSet):
    queryset = FastqRead.objects.all()
    serializer_class = FastqReadSerializer


router = routers.DefaultRouter()
router.register(r'runs', MinionRunViewSet)
router.register(r'reads', FastqReadViewSet)
router.register(r'statistics', RunStatisticViewSet)


urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^api/v1/', include(router.urls)),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
]
