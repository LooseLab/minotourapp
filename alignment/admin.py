from django.contrib import admin

# Register your models here.
from .models import PafRoughCov, PafStore, PafSummaryCov, SamStore


class SamStoreAdmin(admin.ModelAdmin):
    list_display = ('run_id',)
    list_filter = ('run_id',)
    exclude = ('read_id',)
    ordering = ('run_id', )


admin.site.register(SamStore, SamStoreAdmin)
admin.site.register(PafStore)
admin.site.register(PafRoughCov)
admin.site.register(PafSummaryCov)
