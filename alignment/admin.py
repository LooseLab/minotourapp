from django.contrib import admin

from .models import PafRoughCov, PafStore, PafSummaryCov


admin.site.register(PafStore)
admin.site.register(PafRoughCov)
admin.site.register(PafSummaryCov)
