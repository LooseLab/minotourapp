from django.conf.urls import url

from jobs import views as jobs_views
from reads import views


urlpatterns = [
    url(
        r'^api/v1/runs/$',
        views.run_list,
        name="run-list"),
    url(
        r'^api/v1/currentruns/$',
        views.current_run_list,
        name="current-run-list"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9a-z-]+)/$',
        views.run_detail,
        name="run-detail"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/reads/$',
        views.read_list,
        name="fastqread-list"
    ),
    url(
        r'^api/v1/read/$',
        views.read_list_new,
        name="fastqread-list-new"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/readnames/$',
        views.readname_list,
        name="fastqreadnames-list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/yield/$',
        views.cumulative_read_count,
        name="cumulative_read_count"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/rundetails/$',
        views.minION_run_status_list,
        name="minIONrunstatus_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/$',
        views.minION_run_stats_list,
        name="minIONrunstats_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/runstats/(?P<checkid>[0-9]+)/$',
        views.minION_run_stats_latest,
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
        r'^api/v1/runs/(?P<pk>[0-9]+)/histogramsummary/$',
        views.run_histogram_summary,
        name="runhistogramsummary-detail"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/channelsummary/$',
        views.run_channel_summary,
        name="runchannelsummary-list"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/barcodes/$',
        views.barcode_list,
        name="runbarcode-list"),
    url(
        r'^api/v1/reads/(?P<pk>[0-9]+)/$',
        views.read_detail,
        name="fastqread-detail"),
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
        r'^api/v1/activeminions/$',
        views.activeminion_list,
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
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/messagessince/(?P<starttime>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z)/(?P<endtime>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d{3}Z)/$',
        views.sinceminion_messages_list,
        name="sincemessages-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/status/$',
        views.minION_status_list,
        name="minIONstatus-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/$',
        views.minION_scripts_list,
        name="minIONscripts-list"),
    url(
        r'^api/v1/minions/(?P<pk>[0-9]+)/scripts/(?P<nk>[0-9]+)/$',
        views.minION_scripts_detail,
        name="minIONscripts-detail"),
    url(
        r'^api/v1/events/$',
        views.events_type_list,
        name="eventstype-list"),
    url(
        r'^api/v1/events/(?P<pk>[0-9]+)/$',
        views.events_type_detail,
        name="minioneventtype-detail"),
    url(
        r'^api/v1/barcodes/(?P<pk>[0-9]+)/$',
        views.barcode_detail,
        name="barcode-detail"),
    url(
        r'^api/v1/barcodegroupss/(?P<pk>[0-9]+)/$',
        views.barcodegroup_detail,
        name="barcodegroup-detail"),
    url(
        r'^api/v1/tasks/$',
        views.tasks_detail,
        name="tasks-detail"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/tasks/$',
        views.tasks_detail_all,
        name="tasks-detail-all"),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/settask/$',
        views.set_task_detail_all,
        name="set-task-detail-all"),
    url(
        r'^api/v1/flowcells/$',
        views.flowcell_list,
        name="flowcell-list"),
    url(
        r'^api/v1/activeflowcells/$',
        views.flowcell_list_active,
        name="flowcell-list-active"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9A-Za-z-_]+)/$',
        views.flowcell_detail,
        name="flowcell-detail"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcode/$',
        views.flowcell_summary_barcode,
        name="flowcellsummarybarcode-detail"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcodebyminute/$',
        views.flowcell_summary_barcode_by_minute,
        name="flowcellsummarybarcodebyminute-detail"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcodebyminute_bases/$',
        views.flowcell_summary_barcode_by_minute_bases,
        name="flowcellsummarybarcodebyminute-detail_bases"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcodebyminute_quality/$',
        views.flowcell_summary_barcode_by_minute_quality,
        name="flowcellsummarybarcodebyminute-detail_quality"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcodebyminute_length/$',
        views.flowcell_summary_barcode_by_minute_length,
        name="flowcellsummarybarcodebyminute-detail_length"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcodebyminute_maxlength/$',
        views.flowcell_summary_barcode_by_minute_maxlength,
        name="flowcellsummarybarcodebyminute-detail_maxlength"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcodebyminute_speed/$',
        views.flowcell_summary_barcode_by_minute_speed,
        name="flowcellsummarybarcodebyminute-detail_speed"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/histogramsummary/$',
        views.flowcell_histogram_summary,
        name="flowcellhistogramsummary-detail"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/channelsummary/$',
        views.flowcell_channel_summary,
        name="flowcellchannelsummary-list"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/channelsummary_readcount/$',
        views.flowcell_channel_summary_readcount,
        name="flowcellchannelsummary-list-readcount"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/channelsummary_readkb/$',
        views.flowcell_channel_summary_readkb,
        name="flowcellchannelsummary-list-readkb"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/rundetails/$',
        views.flowcell_run_status_list,
        name="flowcellrunstatus_list"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/runstats/$',
        views.flowcell_run_stats_list,
        name="flowcellrunstats_list"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/runstats/(?P<checkid>[0-9]+)/$',
        views.flowcell_run_stats_latest,
        name="flowcellrunstats_latest"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/tasks/$',
        views.flowcell_tasks_detail_all,
        name="flowcelltasks-detail-all"),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/settask/$',
        views.flowcellset_task_detail_all,
        name="flowcellset-task-detail-all"),
    url(
        r'^api/v1/tabs/(?P<pk>[0-9]+)/$',
        views.tabs_details,
        name="tabs-details"
    ),
    url(
        r'^api/v1/flowcells/tabs/(?P<pk>[0-9]+)/$',
        views.flowcell_tabs_details,
        name="tabs-details"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/minknow-messages/$',
        views.minknow_message_list_by_flowcell,
        name="minknow-message-list-by-flowcell"
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
    url(
        r'^api/v1/jobtype/$',
        jobs_views.job_type_list,
        name='jobtype-list'
    ),
    url(
        r'^api/v1/jobtype/(?P<pk>[0-9A-Za-z-]+)/$',
        jobs_views.job_type_detail,
        name='jobtype-detail'
    ),
]
