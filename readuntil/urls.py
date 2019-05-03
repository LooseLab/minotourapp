from django.conf.urls import url

from readuntil import views

urlpatterns = [
    url(
        r'^api/v1/readuntil/benefitdata/master$',
        views.get_benefit_data,
        name="benefit-data"
    ),
    url(
        r'^api/v1/readuntil/benefitdata/detail/$',
        views.get_benefit_data_complete,
        name="benefit-data-non-smoothed"
    ),
]
