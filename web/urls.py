from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name="run_list"),
    url(r'^current/$', views.current, name="current_run"),
    url(r'^log_in/$', views.log_in, name="log_in"),
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^private/index/$', views.private_index, name='private-index'),
    url(r'^private/minup/$', views.minup, name='minup'),
    url(r'^private/profile/$', views.profile, name='profile'),
    url(r'^private/tutorial/$', views.tutorial, name='tutorial'),
    url(r'^private/previous-runs/$', views.previous_run, name='previous-runs'),
    url(r'^private/previous-flowcells/$', views.previous_flowcell, name='previous-flowcells'),
    url(r'^private/current-runs/$', views.current_run, name='current-runs'),
    url(r'^private/runs/(?P<pk>[0-9]+)/$', views.run_index, name='run-index'),
    url(r'^private/flowcells/(?P<pk>[0-9]+)/$', views.flowcell_index, name='flowcell-index'),
    url(r'^private/remote-control/$', views.remotecontrol, name='remotecontrol'),
    url(r'^private/sandbox/$', views.sandbox, name='sandbox'),
]
