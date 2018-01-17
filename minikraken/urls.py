from django.conf.urls import url

from minikraken import views

urlpatterns = [
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/krakenraw/$',
        views.kraken_raw_list,
        name="kraken_raw_list"
        #(?P<bc>[0-9]+)/(?P<ch>[0-9]+)/
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/krakenparse/$',
        views.kraken_parse_list,
        name="kraken_parse_list"
        #(?P<bc>[0-9]+)/(?P<ch>[0-9]+)/
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/krakenraw/$',
        views.flowcellkraken_raw_list,
        name="flowcellkraken_raw_list"
        #(?P<bc>[0-9]+)/(?P<ch>[0-9]+)/
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/krakenparse/$',
        views.flowcellkraken_parse_list,
        name="flowcellkraken_parse_list"
        #(?P<bc>[0-9]+)/(?P<ch>[0-9]+)/
    ),
    ]