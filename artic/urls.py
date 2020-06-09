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
    )
]