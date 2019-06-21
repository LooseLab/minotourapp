from django.conf.urls import url

from readuntil import views

urlpatterns = [
    url(
        r'^api/v1/readuntil/benefitdata/master$',
        views.get_benefit_data_master,
        name="benefit-data"
    ),
    url(
        r'^api/v1/readuntil/benefitdata/detail/$',
        views.get_benefit_data_detail,
        name="benefit-data-detail"
    ),
    url(
        r'^api/v1/readuntil/benefitdata/chromosomes/$',
        views.get_chromosomes,
        name="benefit-data-chromosomes"
    )
]
