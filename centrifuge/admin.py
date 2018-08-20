from django.contrib import admin

from centrifuge.models import CentOutput, LineageValues, MetaGenomicsMeta, SankeyLinks

admin.site.register(CentOutput)
admin.site.register(LineageValues)
admin.site.register(MetaGenomicsMeta)
admin.site.register(SankeyLinks)
