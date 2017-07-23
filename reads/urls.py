from django.conf.urls import url

from reads import views

urlpatterns = [
    url(r'^api/v1/runs/$', views.run_list, name="run-list"),
    url(r'^api/v1/currentruns/$', views.current_run_list, name="current-run-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/$', views.run_detail, name="minionrun-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/reads/$', views.read_list, name="fastqread-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/readnames/$', views.readname_list, name="fastqreadnames-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/yield/$', views.cumulative_read_count, name="cumulative_read_count"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/rundetails/$', views.minION_run_status_list, name="minIONrunstatus_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/$', views.minION_run_stats_list, name="minIONrunstats_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/(?P<checkid>[0-9]+)/$', views.minION_run_stats_latest, name="minIONrunstats_latest"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/summary/$', views.run_summary, name="runsummary-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/summarybyminute/$', views.run_summary_by_minute, name="runsummarybyminute-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/summarybyminute/(?P<last>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)/$', views.run_summary_by_minute, name="runsummarybyminute-detail2"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/summarybarcode/$', views.run_summary_barcode, name="runsummarybarcode-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/summarybarcodebyminute/$', views.run_summary_barcode_by_minute, name="runsummarybarcodebyminute-detail"),
    url(r'^api/v1/reads/(?P<pk>[0-9]+)/$', views.read_detail, name="fastqread-detail"),
    url(r'^api/v1/readtypes/$', views.read_type_list, name="fastqreadtype-list"),
    url(r'^api/v1/readtypes/(?P<pk>[0-9]+)/$', views.read_type_detail, name="fastqreadtype-detail"),
    url(r'^api/v1/minions/$', views.minion_list, name="minION_list"),
    url(r'^api/v1/activeminions/$', views.activeminion_list, name="activeminION_list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/$', views.minION_detail, name="minion-detail"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/events/$', views.minION_events_list, name="minIONevents-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/control/$', views.minION_control_list, name="minIONcontrol-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/control/(?P<checkid>[0-9]+)/$', views.minION_control_update, name="minIONcontrol-update"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/messages/$', views.minION_messages_list, name="minIONmessages-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/recentmessages/$', views.recentminION_messages_list, name="recentmessages-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/status/$', views.minION_status_list, name="minIONstatus-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/$', views.minION_scripts_list, name="minIONscripts-list"),
    url(r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/(?P<nk>[0-9]+)/$', views.minION_scripts_detail, name="minIONscripts-detail"),
    url(r'^api/v1/events/$', views.events_type_list, name="eventstype-list"),
    url(r'^api/v1/events/(?P<pk>[0-9]+)/$', views.events_type_detail, name="minioneventtype-detail"),
]
