from django.conf.urls import url

from reads import views

urlpatterns = [
    url(r'^api/v1/runs/$', views.run_list, name="run_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/$', views.run_detail, name="minionrun-detail"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/reads/$', views.read_list, name="fastqread-list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/readnames/$', views.readname_list, name="fastqreadnames-list"),
    url(r'^api/v1/reads/(?P<pk>[0-9]+)/$', views.read_detail, name="fastqread-detail"),
    url(r'^api/v1/readtypes/$', views.read_type_list, name="fastqreadtype-list"),
    url(r'^api/v1/readtypes/(?P<pk>[0-9]+)/$', views.read_type_detail, name="fastqreadtype-detail"),
]
