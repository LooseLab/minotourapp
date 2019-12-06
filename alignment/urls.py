from django.conf.urls import url

from alignment import views

urlpatterns = [
    # url(
    #     r'^api/v1/flowcells/(?P<flowcell_id>[0-9]+)/pafcover/(?P<barcodegroup_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/(?P<start>-?[0-9]+)/(?P<end>[0-9]+)/$',
    #     views.flowcell_paf_alignment_list,
    #     name="flowcell_paf_alignment_list"
    # ),
    # url(
    #     r'^api/v1/flowcells/(?P<flowcell_id>[0-9]+)/pafcover/(?P<barcodegroup_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/$',
    #     views.rough_coverage_complete_chromosome_flowcell,
    #     name="rough_coverage_complete_chromosome_flowcell"
    # ),
    url(
        r'^api/v1/pafcoverage/(?P<task_id>[0-9]+)/(?P<barcode_name>[a-zA-Z0-9 ]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/$',
        views.rough_coverage_complete_chromosome_flowcell,
        name="rough_coverage_complete_chromosome_flowcell"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/pafsummary/$',
        views.flowcell_paf_summary_cov,
        name="flowcell_paf_summary_cov"
    ),
    url(
        r'^api/v1/flowcells/(?P<flowcell_id>[0-9]+)/references/$',
        views.flowcellreferences_used_by_run,
        name="flowcellreferences_used_by_run"
    ),

    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/pafsummarytable/$',
        views.paf_summary_json,
        name="flowcell_paf_summary"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/chromosomecoverage/$',
        views.get_coverage_summary,
        name="get_coverage_summary"
    )
]
