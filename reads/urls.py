from django.conf.urls import url

from reads import views

urlpatterns = [
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/barcode-proportion/$',
        views.proportion_of_total_reads_in_barcode_list,
        name="barcode-proportion"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/minknow-messages-html/$',
        views.minknow_message_html,
        name="minknow-messages-html"
    ),
    url(
        r'^api/v1/flowcells/reactivate/$',
        views.reactivate_flowcell,
        name="reactivate_flowcell"
    ),
    url(
        r'^api/v1/runs/$',
        views.run_list,
        name="run-list"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9a-z-]+)/$',
        views.run_detail,
        name="run-detail"
    ),
    # url(
    #     r'^api/v1/runs/(?P<pk>[0-9]+)/reads/$',
    #     views.read_list,
    #     name="fastqread-list"
    # ),
    url(
        r'^api/v1/read/$',
        views.read_list_new,
        name="fastqread-list-new"
    ),
    url(
        r'^api/v1/readsextra/$',
        views.readextra_list,
        name="readextra_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/readnames/$',
        views.readname_list,
        name="fastqreadnames-list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/rundetails/$',
        views.list_minion_run_status,
        name="minIONrunstatus_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/$',
        views.minion_run_stats_list,
        name="minIONrunstats_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/(?P<checkid>[0-9]+)/$',
        views.minion_run_stats_latest,
        name="minIONrunstats_latest"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/summarybarcode/$',
        views.run_summary_barcode,
        name="runsummarybarcode-detail"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/summarybarcodebyminute/$',
        views.run_summary_barcode_by_minute,
        name="runsummarybarcodebyminute-detail"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9A-Za-z-_]+)/files/$',
        views.fastq_file,
        name="fastqfile-list",
    ),
    url(
        # r'^api/v1/minions/(?P<pk>[0-9A-Za-z-_]+)/$',
        r'^api/v1/files/(?P<pk>[0-9]+)/$',
        views.fastq_detail,
        name="fastqfile-detail"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/barcodes/$',
        views.barcode_list,
        name="runbarcode-list"),
    url(
        r'^api/v1/readtypes/$',
        views.read_type_list,
        name="fastqreadtype-list"),
    url(
        r'^api/v1/readtypes/(?P<pk>[0-9]+)/$',
        views.read_type_detail,
        name="fastqreadtype-detail"),
    url(
        r'^api/v1/minions/$',
        views.minion_list,
        name="minION_list"),
    url(
        r'^api/v1/active_minions/$',
        views.active_minion_list,
        name="activeminION_list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9A-Za-z-_]+)/$',
        views.minION_detail,
        name="minion-detail"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/events/$',
        views.minION_events_list,
        name="minIONevents-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/control/$',
        views.minION_control_list,
        name="minIONcontrol-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/control/(?P<checkid>[0-9]+)/$',
        views.minION_control_update,
        name="minIONcontrol-update"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/messages/$',
        views.minion_messages_list,
        name="minIONmessages-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/recentmessages/$',
        views.recentminion_messages_list,
        name="recentmessages-list"),
    # url( TODO unused
    #     r'^api/v1/minions/(?P<pk>[0-9]+)/messagessince/(?P<starttime>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)/(?P<endtime>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d{3}Z)/$',
    #     views.sinceminion_messages_list,
    #     name="sincemessages-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/status/$',
        views.list_minion_info,
        name="minion-info-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/$',
        views.minION_scripts_list,
        name="minIONscripts-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/(?P<nk>[0-9]+)/$',
        views.minION_scripts_detail,
        name="minIONscripts-detail"),
    url(
        r'^api/v1/minioneventtypes/$',
        views.events_type_list,
        name="eventstype-list"),
    url(
        r'^api/v1/minioneventtypes/(?P<pk>[0-9]+)/$',
        views.events_type_detail,
        name="minioneventtype-detail"),
    url(
        r'^api/v1/barcodes/(?P<pk>[0-9]+)/$',
        views.barcode_detail,
        name="barcode-detail"),
    url(
        r'^api/v1/flowcells/$',
        views.flowcell_list,
        name="flowcell-list"),
    url(
        r'^api/v1/flowcells/barcodes/$',
        views.flowcell_barcodes_list,
        name="flowcell-barcodes"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9A-Za-z-_]+)/$',
        views.flowcell_detail,
        name="flowcell-detail"),
    url(  # TODO delete?
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcode/$',
        views.flowcell_summary_barcode,
        name="flowcellsummarybarcode-detail"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/statistics/$',
        views.flowcell_statistics,
        name="flowcell-statistics"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/histogram-summary/$',
        views.flowcell_histogram_summary,
        name="flowcell-histogram-summary-detail"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/channel-summary/$',
        views.flowcell_channel_summary,
        name="flowcell-channel-summary"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/runstats/(?P<checkid>[0-9]+)/$',
        views.flowcell_run_stats_latest,
        name="flowcellrunstats_latest"
    ),
    # TODO IS THIS VIEW USED
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/tasks/$',
        views.flowcell_tasks_detail_all,
        name="flowcelltasks-detail-all"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/sharing/$',
        views.flowcell_sharing,
        name="flowcell-sharing"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/sharing/delete/$',
        views.flowcell_sharing_delete,
        name="flowcell-sharing-delete"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/tabs/$',
        views.flowcell_tabs_list,
        name="tabs-details"
    ),
    url(
        r'^api/v1/grouprun/$',
        views.grouprun_list,
        name="grouprun-list"
    ),
    url(
        r'^api/v1/grouprun/(?P<pk>[0-9A-Za-z-_]+)/$',
        views.grouprun_detail,
        name="grouprun-detail"
    ),
    url(
        r'^api/v1/grouprun-membership/$',
        views.grouprun_membership_list,
        name="grouprun-membership-list"
    ),
    url(
        r'^api/v1/barcode/$',
        views.barcode_list_new,
        name="barcode-list-new"
    ),
    url(
        r'^api/v1/version/$',
        views.version,
        name='version'
    ),
    #
    # The links below generated html
    #
    url(
        r'^flowcells/(?P<pk>[0-9]+)/run_summaries_html/$',
        views.flowcell_run_summaries_html,
        name="flowcell_run_summaries_html"
    ),

    url(
        r'^flowcells/(?P<pk>[0-9]+)/flowcell_run_basecalled_summary_html/$',
        views.flowcell_run_basecalled_summary_html,
        name="flowcell_run_basecalled_summary_html"
    ),
    url(
        r'^flowcells/(?P<pk>[0-9]+)/basecalled-summary-html/$',
        views.flowcell_basecalled_summary_html,
        name="flowcell-summary-html"
    ),
    url(
        r'^api/v1/tasks/action$',
        views.task_control,
        name="restart-task",
    ),

    url(
        r'^api/v1/tasktypes/$',
        views.task_types_list,
        name="task-types-list"),
    url(
        r'^api/v1/tasks/$',
        views.job_master_list,
        name='task_list'
    ),
    # TODO web only
    url(
        r'^api/v1/active_minions/run_stats/$',
        views.remote_control_run_stats,
        name="remote-control-run_stats"
    ),

]
