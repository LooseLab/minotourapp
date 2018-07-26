from django.conf.urls import url

from readuntil import views

urlpatterns = [
    url(
        r'^api/v1/rejectedfastread/$',
        views.rejectedfastqread_list,
        name="rejectedfastqread-list"),
]
