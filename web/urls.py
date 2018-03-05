from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name="run_list"),
    url(r'^log_in/$', views.log_in, name="log_in"),
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^private/index/$', views.private_index, name='private-index'),
    url(r'^private/minup/$', views.minup, name='minup'),
    url(r'^private/profile/$', views.profile, name='profile'),
    url(r'^private/tutorial/$', views.tutorial, name='tutorial'),
    url(r'^private/runs/$', views.runs, name='runs'),
    url(r'^private/flowcells/$', views.flowcells, name='flowcells'),
    url(r'^private/runs/(?P<pk>[0-9]+)/$', views.run_index, name='run-index'),
    url(r'^private/flowcells/(?P<pk>[0-9]+)/$', views.flowcell_index, name='flowcell-index'),
    url(r'^private/remote-control/$', views.remotecontrol, name='remotecontrol'),
    url(r'^private/sandbox/$', views.sandbox, name='sandbox'),
]
