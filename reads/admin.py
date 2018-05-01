from django.contrib import admin

from .models import (Barcode, FastqRead, FastqReadType,
                     MinIONControl, MinIONEvent, MinIONEventType,
                     MinIONmessages, MinIONRunStats, MinIONRunStatus,
                     MinIONScripts, MinIONStatus, Run, UserOptions, ChannelSummary, HistogramSummary,
                     RunStatisticBarcode, RunSummaryBarcode, GroupRun, FlowCellRun, GroupBarcode)


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
admin.site.register(ChannelSummary)
admin.site.register(HistogramSummary)
admin.site.register(UserOptions)
admin.site.register(Barcode)
admin.site.register(GroupBarcode)
admin.site.register(FlowCellRun)
admin.site.register(GroupRun)
