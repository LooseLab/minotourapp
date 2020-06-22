from django.conf.urls import url

from communication import views

urlpatterns = [
    url(
        r'^api/v1/communication/messages/conditions$',
        views.get_create_delete_conditions,
        name="conditions"
    ), url(
        r'^api/v1/communication/(?P<pk>[0-9]+)/chromosomecoverage/$',
        views.get_coverage_summary,
        name="get_coverage_summary"
    )
]
