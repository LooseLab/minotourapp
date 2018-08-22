from django.conf.urls import url
from rest_framework.routers import DefaultRouter

from centrifuge.views import CentViewSet, start_centrifuge_view, DefaultViewSet, \
    get_or_set_cartmap, cent_sankey_two, vis_table_data, metaview


router = DefaultRouter()
router.register(prefix="cent", viewset=CentViewSet, base_name="cent")
router.register(prefix="defaulttargets", viewset=DefaultViewSet, base_name="default")
# TODO Remove the three lines above.

urlpatterns = [
    url(r"^start/$", start_centrifuge_view),
    url("sankey", cent_sankey_two),
    url("cartmaptargets", get_or_set_cartmap),
    url("donut", vis_table_data),
    url("table", vis_table_data),
    url("metaview", metaview)
]
