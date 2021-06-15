from django.conf.urls import url

from alignment import views

urlpatterns = [
    url(
        r'^api/v1/alignment/coverage/(?P<task_id>[0-9]+)/(?P<barcode_id>[a-zA-Z0-9 ]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/(?P<rejection_id>[0-9]+)/$',
        views.rough_coverage_complete_chromosome_flowcell,
        name="rough_coverage_complete_chromosome_flowcell"
    ),
    url(
        r'^api/v1/alignment/coverage/(?P<task_id>[0-9]+)/(?P<barcode_id>[a-zA-Z0-9 ]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/(?P<rejection_id>[0-9]+)/(?P<start>[0-9]+)/(?P<end>[0-9]+)/$',
        views.rough_coverage_partial_chromosome_flowcell,
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
    )
]
