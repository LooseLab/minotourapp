from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name="run_list"),
    url(r'^current/$', views.current, name="current_run"),
    url(r'^log_in/$', views.log_in, name="log_in"),
    url(r'^private/index/$', views.private_index, name='private-index'),
    url(r'^private/external-links/$', views.external_links, name='external-links'),
    url(r'^private/minup/$', views.minup, name='minup'),
    url(r'^private/profile/$', views.profile, name='profile'),
    url(r'^private/tutorial/$', views.tutorial, name='tutorial'),
    url(r'^private/previous-runs/$', views.previous_run, name='previous-runs'),
    url(r'^private/current-runs/$', views.current_run, name='current-runs'),
    url(r'^private/runs/(?P<pk>[0-9]+)/$', views.run_index, name='run-index'),
    url(r'^private/remote-control/$', views.remotecontrol, name='remotecontrol'),
    url(r'^private/prev-remote-control/$', views.prevremotecontrol, name='prevremotecontrol'),
]
