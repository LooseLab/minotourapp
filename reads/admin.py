from django.contrib import admin

from jobs.models import JobMaster, JobType
from stats.models import (ChannelSummary, HistogramSummary,
                          RunStatisticBarcode, RunSummaryBarcode)

from .models import (Barcode, BarcodeGroup, FastqRead, FastqReadType,
                     FlowCellRun, MinIONControl, MinIONEvent, MinIONEventType,
                     MinIONmessages, MinIONRunStats, MinIONRunStatus,
                     MinIONScripts, MinIONStatus, Run, UserOptions)


class RunStatisticsBarcodeAdmin(admin.ModelAdmin):
    list_display = ('run_id', 'type')
    list_filter = ('run_id', )
    ordering = ('run_id', 'type__name')


class RunSummaryBarcodeAdmin(admin.ModelAdmin):
    list_display = ('run_id', 'barcode', 'type')
    list_filter = ('run_id', )
    ordering = ('run_id', 'barcode', 'type')


admin.site.register(FastqRead)
admin.site.register(FastqReadType)
admin.site.register(Run)
admin.site.register(MinIONEventType)
admin.site.register(MinIONEvent)
admin.site.register(MinIONScripts)
admin.site.register(MinIONStatus)
admin.site.register(MinIONRunStatus)
admin.site.register(MinIONRunStats)
admin.site.register(MinIONmessages)
admin.site.register(MinIONControl)
admin.site.register(RunStatisticBarcode, RunStatisticsBarcodeAdmin)
admin.site.register(RunSummaryBarcode, RunSummaryBarcodeAdmin)
admin.site.register(JobMaster)
admin.site.register(JobType)
admin.site.register(ChannelSummary)
admin.site.register(HistogramSummary)
admin.site.register(UserOptions)
admin.site.register(Barcode)
admin.site.register(BarcodeGroup)
admin.site.register(FlowCellRun)
