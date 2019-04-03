from django.conf.urls import url

from readuntil import views

urlpatterns = [
    url(
        r'^api/v1/rejectedfastread/$',
        views.rejectedfastqread_list,
        name="rejectedfastqread-list"),
    url(
        r'^api/v1/readuntil/benefitdata/$',
        views.get_benefit_data,
        name="benefit-data"
    ),
    url(
        r'^api/v1/readuntil/benefitdata/complete/$',
        views.get_benefit_data_complete,
        name="benefit-data-non-smoothed"
    ),
]
