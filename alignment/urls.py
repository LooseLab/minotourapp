from django.conf.urls import url

from alignment import views

urlpatterns = [
    url(
        r'^api/v1/runs/(?P<run_id>[0-9]+)/pafcover/(?P<barcode_id>[0-9]+)/(?P<read_type_id>[0-9]+)/(?P<chromosome_id>[0-9]+)/$',
        views.paf_alignment_list,
        name="paf_alignment_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafsummary/$',
        views.paf_alignment_summary,
        name="paf_alignment_summary"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafsummarytrans/$',
        views.paf_alignment_transcriptome_summary,
        name="paf_alignment_transcriptome_summary"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/(?P<ch>[0-9]+)/$',
        views.paf_alignment_ch_list,
        name="paf_alignment_ch_list"
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/barcode/(?P<bc>[0-9]+)/chromosome/(?P<ch>[0-9]+)/type/(?P<ty>[0-9]+)/pos/(?P<po>[0-9]+)/len/(?P<ln>[0-9]+)/$',
        views.paf_test,
        name="paf_test"
    ),
    #url(
    #    r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/chromosome/$',
    #)
]
