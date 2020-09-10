from django.contrib import admin

from .models import PafRoughCov, PafSummaryCov, PafRoughCovIntermediate

admin.site.register(PafRoughCov)
admin.site.register(PafSummaryCov)
admin.site.register(PafRoughCovIntermediate)
