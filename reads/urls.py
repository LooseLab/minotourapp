from django.conf.urls import url

from reads import views

urlpatterns = [
    #url(r'^api/v1/dataset/$', views.data_set, name="data_set"),
    url(r'^api/v1/runs/$', views.run_list, name="run_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/$', views.run_detail, name="minionrun-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/reads/$', views.read_list, name="fastqread-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/readnames/$', views.readname_list, name="fastqreadnames-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/yield/$', views.cumulative_read_count, name="cumulative_read_count"),
    url(r'^api/v1/reads/(?P<pk>[0-9]+)/$', views.read_detail, name="fastqread-detail"),
    url(r'^api/v1/readtypes/$', views.read_type_list, name="fastqreadtype-list"),
    url(r'^api/v1/readtypes/(?P<pk>[0-9]+)/$', views.read_type_detail, name="fastqreadtype-detail"),
    url(r'^api/v1/minions/$', views.minion_list, name="minION_set"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/$', views.minION_event, name="minIONevent"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/events/$', views.minION_events_list, name="minIONevents-list"),
    url(r'^api/v1/minions/events/$', views.minION_events_type, name="minIONeventstype-list")
]
