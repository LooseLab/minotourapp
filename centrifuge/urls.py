from django.conf.urls import url

from centrifuge import views

from centrifuge.views import get_or_set_cartmap, cent_sankey_two, vis_table_or_donut_data, metaview, get_target_mapping

urlpatterns = [
    # Return all the data for the sankey diagram in the correct form, called from sankey.js in the web app
    url(r"^sankey/$", cent_sankey_two),
    # TODO currently unused, might become part of the mapping step
    url(r"^cartmaptargets/$", get_or_set_cartmap),
    # Return the data for the donut chart, in the correct form, called form donutdothat.js
    url(r"^donut/$", vis_table_or_donut_data),
    # Return the data for the total reads table, called from AllResultsTable.js
    url(r"^table/$", vis_table_or_donut_data),
    # Return the data for the header containing metadata at the top of the vis page
    url(r"^metaview/$", metaview),
    # Return the data on for the target table
    url(r"^mapped_targets/$", get_target_mapping),
    url(r'^api/v1/flowcells/(?P<pk>[0-9]+)/metagenomic_barcodes/$', views.metagenomic_barcodes),
]
