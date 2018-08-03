from django.conf.urls import url
from rest_framework.routers import DefaultRouter
from centrifuge.views import CentViewSet, cent_donut, start_centrifuge_view, DefaultViewSet, \
    get_or_set_cartmap, RefViewSet, cent_sankey, vis_table_data, metaview

# django rest framework router
# router = DefaultRouter()
# router.register(prefix="cent", viewset=CentViewSet, base_name="cent")
# router.register(prefix="defaulttargets", viewset=DefaultViewSet, base_name="default")
# router.register(prefix="ref", viewset=RefViewSet, base_name="references")

urlpatterns = [
    url(r"^start/$", start_centrifuge_view),
    url(r"^sankey", cent_sankey),
    url(r"^cartmaptargets", get_or_set_cartmap),
    url(r"^donut", cent_donut),
    url(r"^table", vis_table_data),
    url(r"^metaview", metaview)
]

# urlpatterns += router.urls
