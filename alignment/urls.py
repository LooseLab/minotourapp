from django.conf.urls import url

from alignment import views

urlpatterns = [
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/$',
        views.paf_alignment_list,
        name="paf_alignment_list"
    ),
    #url(
    #    r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/chromosome/$',
    #)
]