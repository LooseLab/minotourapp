"""
Admin.py
"""
from django.contrib import admin

from metagenomics.models import CentrifugeOutput, LineageValue, Metadata, SankeyLink, \
    MappingResult, TargetMappedReadId, MappingTarget

# Register all the Models to the django admin site for this app
admin.site.register(CentrifugeOutput)
admin.site.register(LineageValue)
admin.site.register(Metadata)
admin.site.register(SankeyLink)
admin.site.register(MappingResult)
admin.site.register(TargetMappedReadId)
admin.site.register(MappingTarget)


