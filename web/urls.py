from django.conf.urls import url

from . import views
urlpatterns = [
    url(r'^signup/$', views.signup, name='signup'),
    url(r'^private/index/$', views.private_index, name='private-index'),
    url(r'^private/minup/$', views.minup, name='minup'),
    url(r'^private/profile/$', views.profile, name='profile'),
    url(
        r'^private/message/$',
        views.message_details,
        name="message_details"
    ),
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
    url(r'^private/flowcell_manager/$', views.flowcell_manager, name='flowcell_manager'),
    #
    # The links below generate CSV
    #
    url(
        r'^private/flowcells/(?P<pk>[0-9]+)/flowcell_run_stats_csv/$',
        views.flowcell_run_stats_download,
        name="flowcell_run_stats_download"
    ),
    # Return the netagenomics data in CSV format
    url(r'^private/flowcells/(?P<pk>[0-9]+)/flowcell_metagenomics_csv/$', views.metagenomics_data_download,
        name="metagenomics_data_download")

]
