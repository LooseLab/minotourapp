from django.conf.urls import url

from reference import views

urlpatterns = [
    url(
        r'^api/v1/reference/$',
        views.reference_list,
        name="reference-info"
    ),
    url(
        r'^api/v1/reference/upload/$',
        views.reference_list,
        name="reference_upload"
    ),
    url(
        r'^api/v1/reference/manager/$',
        views.reference_manager,
        name="reference-manager"
    )
]
