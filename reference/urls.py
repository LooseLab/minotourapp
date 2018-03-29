from django.conf.urls import url

from reference import views

urlpatterns = [
    url(
        r'^api/v1/reference/$',
        views.reference_all,
        name="referenceinfo"
    ),
    url(
        r'^api/v1/reference/(?P<pk>[0-9]+)/$',
        views.reference_list,
        name="referenceinfo-detail"
    ),
    url(
        r'^api/v1/chromosome/(?P<pk>[0-9]+)/$',
        views.reference_line,
        name="referenceline-detail"
    ),
    #url(
    #    r'^api/v1/runs/(?P<pk>[0-9]+)/pafcover/chromosome/$',
    #)
]
