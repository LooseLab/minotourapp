from django.conf.urls import url

from alignment import views

urlpatterns = [
    url(
        r'^api/v1/alignment/coverage/(?P<task_id>[0-9]+)/(?P<barcode_id>[a-zA-Z0-9 ]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/$',
        views.coverage_master,
        name="rough_coverage_complete_chromosome_flowcell"
    ),
    url(
        r'^api/v1/alignment/coverage/(?P<task_id>[0-9]+)/(?P<barcode_id>[a-zA-Z0-9 ]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/(?P<start>[0-9]+)/(?P<end>[0-9]+)/$',
        views.coverage_detail,
        name="rough_coverage_partial_chromosome_flowcell"
    ),
    url(
        r'^api/v1/alignment/(?P<pk>[0-9]+)/pafsummary/$',
        views.per_barcode_coverage_summary,
        name="per_barcode_coverage_summary"
    ),
    url(
        r'^api/v1/alignment/(?P<flowcell_id>[0-9]+)/mapped-references/$',
        views.mapped_references_by_flowcell_list,
        name="mapped_references_by_flowcell_list"
    ),

    url(
        r'^api/v1/alignment/(?P<pk>[0-9]+)/paf-summary-table/$',
        views.paf_summary_table_json,
        name="flowcell_paf_summary"
    ),
    url(
        r'^api/v1/alignment/(?P<flowcell_pk>[0-9]+)/genome-coverage-summary',
        views.per_genome_coverage_summary,
        name="genome-coverage-summary"
    ),
    url(
        r'^api/v1/alignment/(?P<pk>[0-9]+)/cnv-chart/(?P<barcode_pk>[0-9]+)/(?P<job_pk>[0-9]+)/(?P<expected_ploidy>[0-9])',
        views.cnv_chart,
        name="cnv-chart"
    ),
    url(
        r'^api/v1/alignment/(?P<pk>[0-9]+)/cnv-barcode',
        views.get_cnv_barcodes,
        name="cnv-barcodes"
    ),
    url(
        r'^api/v1/alignment/(?P<pk>[0-9]+)/cnv-chart-detail/(?P<barcode_pk>[0-9]+)/(?P<job_pk>[0-9]+)/(?P<contig_name>[a-zA-Z0-9 ]+)/(?P<exp_ploidy>[0-9]+)/(?P<pen_value>[0-9]+)/(?P<min_diff>[0-9]+)/(?P<bin_slice>[0-9]+)/(?P<median_bin_value>[0-9]+)$',
        views.cnv_detail_chart,
        name="cnv-detail-chart"
    )
]
