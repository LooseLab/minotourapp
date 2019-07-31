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
    ),
    url(
        r'^api/v1/readuntil/rejected_ids/$',
        views.get_or_create_rejected_reads,
        name="rejected_ids"
    ),
    url(
        r'api/v1/reduntil/accepted_ids/$',
        views.get_or_create_accepted_reads,
        name="accepted_ids"
    )
]
