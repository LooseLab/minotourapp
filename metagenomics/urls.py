"""
urls.py
"""
from django.conf.urls import url

from metagenomics import views

urlpatterns = [
    # Return all the data for the sankey diagram in the correct form, called from sankey.js in the web app
    url(r"^api/v1/metagenomics/sankey/$", views.centrifuge_sankey),
    # Return the data for the donut chart, in the correct form, called form donutdothat.js
    url(r"^api/v1/metagenomics/donut/$", views.donut_data),
    # Return the data for the total reads table, called from allResultsTable.js
    url(r"^api/v1/metagenomics/(?P<pk>[0-9]+)/table/$", views.all_results_table),
    # Return the data for the header containing metadata at the top of the vis page
    url(r"^api/v1/metagenomics/metagenomics-metadata/$", views.centrifuge_metadata),
    # Return the data on for the target table
    url(r"^api/v1/metagenomics/mapped-targets/$", views.get_target_mapping),
    # Return barcodes that we have metagenomics data for
    url(r'^api/v1/metagenomics/(?P<pk>[0-9]+)/metagenomic-barcodes/$', views.metagenomic_barcodes),
    # Return the data for the simple metagenomics table
    url(r'^api/v1/metagenomics/alerts$', views.simple_target_mappings),
    # Return the data for a choice of target sets
    url(r'^api/v1/metagenomics/targetsets$', views.get_target_sets),
    # Return the simple alert table data
    url(r'^api/v1/metagenomics/(?P<pk>[0-9]+)/simple-alert/(?P<barcode_name>[0-9A-Za-z-_]+)/$', views.super_simple_alert_list),

]
