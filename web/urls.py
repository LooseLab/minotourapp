from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^private/index/$', views.private_index, name='private-index'),
    url(r'^private/minup/$', views.minup, name='minup'),
    url(r'^private/profile/$', views.profile, name='profile'),
    url(r'^private/tutorial/$', views.tutorial, name='tutorial'),
    url(r'^private/flowcells/$', views.flowcells, name='flowcells'),
    url(r'^private/flowcells/(?P<pk>[0-9]+)/$', views.flowcell_index, name='flowcell-index'),
    url(r'^private/flowcells/(?P<pk>[0-9]+)/reads/$', views.flowcell_reads, name='flowcell-reads'),
    url(r'^private/flowcells_reads_data/$', views.flowcell_reads_data, name='flowcell-reads-data'),
    url(r'^private/remote-control/$', views.remotecontrol, name='remotecontrol'),
    url(r'^private/experiments/$', views.ExperimentList.as_view(), name='experiment-list'),
    url(r'^private/experiments/(?P<pk>[0-9]+)/$', views.experiments_update, name='experiment-update'),
    url(r'^private/experiments/create/$', views.experiments_create, name='experiment-create'),
    url(r'^private/experiments/(?P<pk>[0-9]+)/delete/$', views.ExperimentDelete.as_view(), name='experiment-delete'),
]
