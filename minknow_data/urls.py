"""
Urls for the minknow data app
"""
from django.conf.urls import url

from minknow_data import views

urlpatterns = [
    url(
        r"^api/v1/minknow/active_minions/run_stats/$",
        views.minknow_control_run_stats,
        name="remote-control-run_stats",
    ),
    url(
        r"^api/v1/minknow/test-connect/$",
        views.test_connection,
        name="testy-test",
    ),
    url(
        r"^api/v1/minknow/minions/event-types/$",
        views.events_type_list,
        name="eventstype-list",
    ),
    url(
        r"^api/v1/minknow/flowcells/(?P<pk>[0-9]+)/minknow-messages/$",
        views.minknow_message,
        name="minknow-messages",
    ),
    url(
        r"^api/v1/minknow/minions/runs/(?P<pk>[0-9]+)/$",
        views.minion_run_status_list,
        name="minIONrunstatus_list",
    ),
    url(
        r"^api/v1/minknow/minions/runs/(?P<pk>[0-9]+)/run-stats/$",
        views.minion_run_stats_list,
        name="minIONrunstats_list",
    ),
    url(
        r"^api/v1/minknow/minions/runs/(?P<pk>[0-9]+)/run-stats/(?P<checkid>[0-9]+)/$",
        views.minion_run_stats_latest,
        name="minIONrunstats_latest",
    ),
    url(r"^api/v1/minknow/minions/$", views.minion_list, name="minION_list"),
    url(
        r"^api/v1/minknow/active_minions/$",
        views.active_minion_list,
        name="activeminION_list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9A-Za-z-_]+)/$",
        views.minION_detail,
        name="minion-detail",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/events/$",
        views.minION_events_list,
        name="minIONevents-list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/control/$",
        views.minion_control_list,
        name="minIONcontrol-list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/control/(?P<checkid>[0-9]+)/$",
        views.minION_control_update,
        name="minIONcontrol-update",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/messages/$",
        views.minion_messages_list,
        name="minIONmessages-list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/recent-messages/$",
        views.recentminion_messages_list,
        name="recentmessages-list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/status/$",
        views.minion_info_list,
        name="minion-info-list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/scripts/$",
        views.minION_scripts_list,
        name="minIONscripts-list",
    ),
    url(
        r"^api/v1/minknow/minions/(?P<pk>[0-9]+)/scripts/(?P<nk>[0-9]+)/$",
        views.minION_scripts_detail,
        name="minIONscripts-detail",
    ),
    url(r"^api/v1/minknow/version/$", views.version, name="version"),
    url(
        r"^api/v1/minknow/minioneventtypes/(?P<pk>[0-9]+)/$",
        views.events_type_detail,
        name="minioneventtype-detail",
    ),
    url(
        r"^api/v1/minknow/flowcells/(?P<pk>[0-9]+)/runstats/(?P<check_id>[0-9]+)/$",
        views.flowcell_minknow_stats_list,
        name="flowcell-run-stats-list",
    ),

]