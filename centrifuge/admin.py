from django.contrib import admin

from centrifuge.models import CentOutput, LineageValues, Metadata, SankeyLinks, \
    CentOutputBarcoded, MappingResultsBarcoded, MappingResults, RedReadIds, MappingTargets

# Register all the Models to the django admin site for this app
admin.site.register(CentOutput)
admin.site.register(LineageValues)
admin.site.register(Metadata)
admin.site.register(SankeyLinks)
admin.site.register(CentOutputBarcoded)
admin.site.register(MappingResultsBarcoded)
admin.site.register(MappingResults)
admin.site.register(RedReadIds)
admin.site.register(MappingTargets)


