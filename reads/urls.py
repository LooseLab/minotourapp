from django.conf.urls import url

from reads import views

urlpatterns = [
    url(r'^api/v1/runs/$', views.RunListView.as_view(), name="runs_list"),
    url(r'^api/v1/runs/(?P<pk>[0-9]+)/$', views.RunDetailView.as_view(), name="run_detail"),
]
