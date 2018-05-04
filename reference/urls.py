from django.conf.urls import url

from reference import views

urlpatterns = [
    url(
        r'^api/v1/reference/$',
        views.reference_list,
        name="referenceinfo"
    ),
    url(
        r'^api/v1/reference/(?P<pk>[0-9]+)/$',
        views.reference_detail,
        name="referenceinfo-detail"
    ),
    url(
        r'^api/v1/chromosome/(?P<pk>[0-9]+)/$',
        views.reference_line,
        name="referenceline-detail"
    ),
]
