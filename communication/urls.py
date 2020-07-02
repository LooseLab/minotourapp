from django.conf.urls import url

from communication import views

urlpatterns = [
    url(
        r'^api/v1/communication/conditions/create-condition/$',
        views.simple_condition_create,
        name="simple-create"
    ), url(
        r'^api/v1/communication/(?P<pk>[0-9]+)/conditions/references$',
        views.get_references_for_condtions,
        name="get-references"
    ),
    url(
        r'^api/v1/communication/messages/conditions$',
        views.notification_conditions_list,
        name="conditions"
    ), url(
        r'^api/v1/communication/(?P<pk>[0-9]+)/conditions/chromosomecoverage/$',
        views.get_coverage_summary,
        name="get_coverage_summary"
    ), url(
        r'^api/v1/communication/(?P<pk>[0-9]+)/conditions/extant-types/$',
        views.get_extant_condition_types,
        name="get_extant_conditions"
    )
]
