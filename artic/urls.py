from django.conf.urls import url
from artic import views

urlpatterns = [
    url(
        r'^api/v1/artic/visualisation/barcodes/$',
        views.get_artic_barcodes,
        name="get-artic-barcodes"
    ),
    url(
        r'^api/v1/artic/visualisation/master/$',
        views.get_artic_barcodes,
        name="get-artic-barcodes"
    )
]