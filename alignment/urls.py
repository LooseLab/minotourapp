from django.conf.urls import url

from alignment import views

urlpatterns = [
    url(
        r'^api/v1/flowcells/(?P<flowcell_id>[0-9]+)/pafcover/(?P<barcodegroup_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/(?P<start>-?[0-9]+)/(?P<end>[0-9]+)/$',
        views.flowcell_paf_alignment_list,
        name="flowcell_paf_alignment_list"
    ),
    url(
        r'^api/v1/flowcells/(?P<flowcell_id>[0-9]+)/pafcover/(?P<barcodegroup_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/$',
        views.rough_coverage_complete_chromosome_flowcell,
        name="rough_coverage_complete_chromosome_flowcell"
    ),
    url(
        r'^api/v1/flowcells/(?P<pk>[0-9]+)/pafsummary/$',
        views.flowcellpaf_alignment_summary,
        name="flowcellpaf_alignment_summary"
    ),
    url(
        r'^api/v1/flowcells/(?P<flowcell_id>[0-9]+)/references/$',
        views.flowcellreferences_used_by_run,
        name="flowcellreferences_used_by_run"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/(?P<ch>[0-9]+)/$',
        views.paf_alignment_ch_list,
        name="paf_alignment_ch_list"
    ),
    url(
        r'^api/v1/runs/(?P<run_id>[0-9]+)/pafcover/(?P<barcode_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/(?P<start>-?[0-9]+)/(?P<end>[0-9]+)/$',
        views.paf_alignment_list,
        name="paf_alignment_list"
    ),
    url(
        r'^api/v1/runs/(?P<run_id>[0-9]+)/pafcover/(?P<barcode_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/$',
        views.rough_coverage_complete_chromosome_run,
        name="rough_coverage_complete_chromosome_run"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafsummary/$',
        views.paf_alignment_summary,
        name="paf_alignment_summary"
    ),
    url(
        r'^api/v1/runs/(?P<run_id>[0-9]+)/references/$',
        views.references_used_by_run,
        name="references_used_by_run"
    ),
]
