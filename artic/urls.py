from django.conf.urls import url

from artic import views

urlpatterns = [
    url(
        r'^api/v1/artic/visualisation/barcodes/$',
        views.get_artic_barcodes,
        name="get-artic-barcodes"
    ),
    url(
        r'^api/v1/artic/visualisation/master/$',
        views.get_artic_master_chart_data,
        name="get-artic-barcodes"
    ),
    url(
        r'^api/v1/artic/visualisation/detail/$',
        views.get_artic_detail_chart_data,
        name="get-artic-detail"
    ),
    url(
        r'^api/v1/artic/visualisation/column-charts/$',
        views.get_artic_column_chart_data,
        name="get-artic-column"
    ),
    url(
        r'^api/v1/artic/visualisation/summary-table-data/$',
        views.get_artic_summary_table_data,
        name="get-summary-table-data"
    ),
    url(
        r'^api/v1/artic/barcode-metadata/$',
        views.get_artic_barcode_metadata_html,
        name="get-barcode-metadata"
    ),
    url(
        r'^api/v1/artic/manual-trigger/$',
        views.manually_create_artic_command_job_master,
        name="manual-trigger"
    ),
    url(
        r'^api/v1/artic/build-results/$',
        views.get_results_package,
        name="create-results-bundle"
    ),
    url(
        r'^api/v1/artic/fetch-png-html/$',
        views.png_html,
        name="png-html"
    ),
    url(
        r'^api/v1/artic/fetch-amplicon-bands/$',
        views.get_amplicon_bands_for_master,
        name="amplicon-bands"
    ),
    url(
        r'^api/v1/artic/rerun-command/$',
        views.rerun_artic_command,
        name="rerun-command"
    ),
    url(
        r'^api/v1/artic/pie-chart-data/$',
        views.get_artic_pie_chart_data,
        name="pie-chart-data"
    ),
    url(
        r'^api/v1/artic/(?P<pk>[0-9]+)/firing-conditions',
        views.fire_conditions_list,
        name="fire-conditions-data"
    ),
    url(
        r'^api/v1/artic/(?P<pk>[0-9]+)/results-modal',
        views.get_results_modal_html,
        name="results-modal-html"
    )
]