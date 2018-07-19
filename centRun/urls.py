"""djangularCent URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url("", views.home, name="home")
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url("", Home.as_view(), name="home")
Including another URLconf
    1. Import the include() function: from django.urls import include, url
    2. Add a URL to urlpatterns:  url("blog/", include("blog.urls"))
"""
from django.conf.urls import url
from rest_framework.routers import DefaultRouter
from centRun.views import CentViewSet, cent_donut, start_centrifuge_view, DefaultViewSet, \
    get_or_set_cartmap, RefViewSet, cent_sankey, vis_table_data

# django rest framework router
router = DefaultRouter()
router.register(prefix="cent", viewset=CentViewSet, base_name="cent")
router.register(prefix="defaulttargets", viewset=DefaultViewSet, base_name="default")
router.register(prefix="ref", viewset=RefViewSet, base_name="references")

urlpatterns = [
    url(r"^start/$", start_centrifuge_view),
    url("sankey", cent_sankey),
    url("cartmaptargets", get_or_set_cartmap),
    url("donut", cent_donut),
    url("table", vis_table_data),
]

urlpatterns += router.urls
