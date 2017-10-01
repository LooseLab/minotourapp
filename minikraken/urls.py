from django.conf.urls import url

from minikraken import views

urlpatterns = [
    url(
        r'^api/v1/runs/(?P<pk>[0-9]+)/krakenraw/$',
        views.kraken_raw_list,
        name="kraken_raw_list"
        #(?P<bc>[0-9]+)/(?P<ch>[0-9]+)/
    ),
    ]