"""Urls for the read fish/ sword fish integration"""

from django.conf.urls import re_path

from readfish import views

urlpatterns = [
    re_path(r"^swordfish/(?P<run_id>\w+)/toml$", views.TomlDetail.as_view()),
    re_path(r"^swordfish/(?P<run_id>\w+)/chopchop", views.GetToTheChopper.as_view()),
    re_path(r"^swordfish/(?P<run_id>\w+)/validate/(?P<type>\w+)", views.ValidateSwordfishStartup.as_view()),
    re_path(r"^swordfish/test-connect", views.ConnectHere.as_view())
]