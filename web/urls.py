from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.index, name="run_list"),
    url(r'^current/$', views.current, name="current_run"),
    url(r'^log_in/$', views.log_in, name="log_in"),
    url(r'^private/index/$', views.private_index, name='private-index'),
    url(r'^private/external-links/$', views.external_links, name='external-links'),
    url(r'^private/minup/$', views.minup, name='minup'),
]
