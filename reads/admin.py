from django.contrib import admin

from .models import (Barcode, FastqRead, FastqReadType,
                     MinIONControl, MinIONEvent, MinIONEventType,
                     MinionMessage, MinIONRunStats, MinIONRunStatus,
                     MinIONScripts, MinIONStatus, Run, UserOptions, ChannelSummary, HistogramSummary,
                     RunStatisticBarcode, RunSummaryBarcode, GroupRun, FlowcellStatisticBarcode, FlowcellSummaryBarcode,
                     FlowcellTab)


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
admin.site.register(FlowcellStatisticBarcode)
admin.site.register(FlowcellSummaryBarcode)
admin.site.register(FlowcellTab)
admin.site.register(Run)
admin.site.register(MinIONEventType)
admin.site.register(MinIONEvent)
admin.site.register(MinIONScripts)
admin.site.register(MinIONStatus)
admin.site.register(MinIONRunStatus)
admin.site.register(MinIONRunStats)
admin.site.register(MinionMessage)
admin.site.register(MinIONControl)
admin.site.register(RunStatisticBarcode, RunStatisticsBarcodeAdmin)
admin.site.register(RunSummaryBarcode, RunSummaryBarcodeAdmin)
admin.site.register(ChannelSummary)
admin.site.register(HistogramSummary)
admin.site.register(UserOptions)
admin.site.register(Barcode)
admin.site.register(GroupRun)
