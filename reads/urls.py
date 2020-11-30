from django.conf.urls import url

from reads import views

urlpatterns = [
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/barcode-proportion/$",
        views.proportion_of_total_reads_in_barcode_list,
        name="barcode-proportion",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/reactivate/$",
        views.reactivate_flowcell,
        name="reactivate_flowcell",
    ),
    url(r"^api/v1/reads/runs/$", views.run_list, name="run-list"),
    url(
        r"^api/v1/reads/runs/(?P<pk>[0-9a-z-]+)/$", views.run_detail, name="run-detail"
    ),
    url(r"^api/v1/reads/read/$", views.read_list, name="fastqread-list-new"),
    url(r"^api/v1/reads/reads-extra/$", views.readextra_list, name="readextra_list"),
    url(
        r"^api/v1/reads/runs/(?P<pk>[0-9]+)/read-names/$",
        views.readname_list,
        name="fastqreadnames-list",
    ),
    #
    # url(
    #     r"^api/v1/reads/runs/(?P<pk>[0-9]+)/summarybarcode/$",
    #     views.run_summary_barcode,
    #     name="runsummarybarcode-detail",
    # ),
    # url(
    #     r"^api/v1/reads/runs/(?P<pk>[0-9]+)/summarybarcodebyminute/$",
    #     views.run_summary_barcode_by_minute,
    #     name="runsummarybarcodebyminute-detail",
    # ),
    url(
        r"^api/v1/reads/runs/(?P<pk>[0-9A-Za-z-_]+)/files/$",
        views.fastq_file,
        name="fastqfile-list",
    ),
    url(
        r"^api/v1/reads/files/(?P<pk>[0-9]+)/$",
        views.fastq_detail,
        name="fastqfile-detail",
    ),
    url(r"^api/v1/reads/read-types/$", views.read_type_list, name="fastqreadtype-list"),
    url(
        r"^api/v1/reads/readtypes/(?P<pk>[0-9]+)/$",
        views.read_type_detail,
        name="fastqreadtype-detail",
    ),

    url(
        r"^api/v1/reads/barcodes/(?P<pk>[0-9]+)/$",
        views.barcode_detail,
        name="barcode-detail",
    ),
    url(r"^api/v1/reads/flowcells/$", views.flowcell_list, name="flowcell-list"),
    url(
        r"^api/v1/reads/flowcells/barcodes/$",
        views.flowcell_barcodes_list,
        name="flowcell-barcodes",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9A-Za-z-_]+)/$",
        views.flowcell_detail,
        name="flowcell-detail",
    ),
    # url(  # TODO delete?
    #     r"^api/v1/flowcells/(?P<pk>[0-9]+)/summarybarcode/$",
    #     views.flowcell_summary_barcode,
    #     name="flowcellsummarybarcode-detail",
    # ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/statistics/$",
        views.flowcell_statistics,
        name="flowcell-statistics",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/histogram-summary/$",
        views.flowcell_histogram_summary,
        name="flowcell-histogram-summary-detail",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/channel-summary/$",
        views.flowcell_channel_summary,
        name="flowcell-channel-summary",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/tasks/$",
        views.flowcell_tasks_detail_all,
        name="flowcelltasks-detail-all",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/sharing/$",
        views.flowcell_sharing,
        name="flowcell-sharing",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/sharing/delete/$",
        views.flowcell_sharing_delete,
        name="flowcell-sharing-delete",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/tabs/$",
        views.flowcell_tabs_list,
        name="tabs-details",
    ),
    url(r"^api/v1/reads/barcode/$", views.barcode_list, name="barcode-list"),
    #
    # The links below generated html
    #
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/run-summaries-html/$",
        views.flowcell_run_summaries_html,
        name="flowcell_run_summaries_html",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/flowcell-run-basecalled-summary-html/$",
        views.flowcell_run_basecalled_summary_html,
        name="flowcell_run_basecalled_summary_html",
    ),
    url(
        r"^api/v1/reads/flowcells/(?P<pk>[0-9]+)/basecalled-summary-html/$",
        views.flowcell_basecalled_summary_html,
        name="flowcell-summary-html",
    ),
    url(r"^api/v1/reads/tasks/action$", views.task_control, name="restart-task",),
    url(r"^api/v1/reads/task-types/$", views.task_types_list, name="task-types-list"),
    url(r"^api/v1/reads/tasks/$", views.job_master_list, name="task_list"),

]
