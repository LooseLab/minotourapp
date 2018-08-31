from django.conf.urls import url

from centrifuge.views import CentViewSet, start_centrifuge_view, DefaultViewSet, \
    get_or_set_cartmap, cent_sankey_two, vis_table_data, metaview



urlpatterns = [
    url(r"^start/$", start_centrifuge_view),
    url("sankey", cent_sankey_two),
    url("cartmaptargets", get_or_set_cartmap),
    url("donut", vis_table_data),
    url("table", vis_table_data),
    url("metaview", metaview)
]
