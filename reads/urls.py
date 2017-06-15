from django.conf.urls import url

from reads import views

urlpatterns = [
    #url(r'^api/v1/dataset/$', views.data_set, name="data_set"),
    url(r'^api/v1/runs/$', views.run_list, name="run-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/$', views.run_detail, name="minionrun-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/reads/$', views.read_list, name="fastqread-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/readnames/$', views.readname_list, name="fastqreadnames-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/yield/$', views.cumulative_read_count, name="cumulative_read_count"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/rundetails/$', views.minION_run_status_list, name="minIONrunstatus_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/$', views.minION_run_stats_list, name="minIONrunstats_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/reads_statistics/$', views.run_read_statistics, name="minion-run-read-statistics"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/reads_statistics2/$', views.run_read_statistics2, name="minion-run-read-statistics2"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/barcodes/$', views.run_barcodes, name="minion-run-barcodes"),
    url(r'^api/v1/reads/(?P<pk>[0-9]+)/$', views.read_detail, name="fastqread-detail"),
    url(r'^api/v1/readtypes/$', views.read_type_list, name="fastqreadtype-list"),
    url(r'^api/v1/readtypes/(?P<pk>[0-9]+)/$', views.read_type_detail, name="fastqreadtype-detail"),
    url(r'^api/v1/minions/$', views.minion_list, name="minION_list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/$', views.minION_detail, name="minion-detail"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/events/$', views.minION_events_list, name="minIONevents-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/messages/$', views.minION_messages_list, name="minIONmessages-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/status/$', views.minION_status_list, name="minIONstatus-list"),
    #url(r'^api/v1/minions/(?P<pk>[0-9]+)/rundetails/$', views.minION_liverun_list, name="minIONliverun-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/$', views.minION_scripts_list, name="minIONscripts-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/(?P<nk>[0-9]+)/$', views.minION_scripts_detail, name="minIONscripts-detail"),
    url(r'^api/v1/events/$', views.events_type_list, name="eventstype-list"),
    url(r'^api/v1/events/(?P<pk>[0-9]+)/$', views.events_type_detail, name="minioneventtype-detail"),
]
