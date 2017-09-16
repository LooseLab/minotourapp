from django.conf.urls import url

from alignment import views

urlpatterns = [
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/$',
        views.paf_alignment_list,
        name="paf_alignment_list"
        #(?P<bc>[0-9]+)/(?P<ch>[0-9]+)/
    ),
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/(?P<ch>[0-9]+)/$',
        views.paf_alignment_ch_list,
        name="paf_alignment_ch_list"
    ),
    #url(
    #    r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/chromosome/$',
    #)
]