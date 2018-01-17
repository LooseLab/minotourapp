from django.conf.urls import url

from assembly import views

urlpatterns = [
    url(
        r'^api/v1/runs/(?P<run_id>[0-9]+)/assembly/$',
        views.gfa_output_list,
        name="gfa_output_list"
    ),
]
