"""Urls for the read fish/ sword fish integration"""

from django.conf.urls import re_path

from readfish import views

urlpatterns = [
    re_path(r"^swordfish/<str:run_id>/toml$", views.TomlDetail.as_view()),
    re_path(r"^swordfish/<run_id>/chopchop", views.GetToTheChopper.as_view()),
    re_path(r"^swordfish/(<run_id>)/validate", views.ValidateSwordfishStartup.as_view())
]