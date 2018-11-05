from django.contrib import admin

from centrifuge.models import CentOutput, LineageValues, MetaGenomicsMeta, SankeyLinks, \
    CentOutputBarcoded, BarcodedCartographyMapped, CartographyMapped, RedReadIds, CartographyGuide

# Register all the Models to the django admin site for this app
admin.site.register(CentOutput)
admin.site.register(LineageValues)
admin.site.register(MetaGenomicsMeta)
admin.site.register(SankeyLinks)
admin.site.register(CentOutputBarcoded)
admin.site.register(BarcodedCartographyMapped)
admin.site.register(CartographyMapped)
admin.site.register(RedReadIds)
admin.site.register(CartographyGuide)


