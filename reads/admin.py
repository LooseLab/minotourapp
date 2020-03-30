"""
Register models with the admin interface
"""
from django.contrib import admin

from .models import (Barcode, FastqFile, FastqRead, FastqReadType,
                     MinionControl, MinionEvent, MinionEventType,
                     MinionMessage, MinionRunStats, MinionRunInfo,
                     MinionScripts, MinionInfo, Run, RunSummary, UserOptions, ChannelSummary, HistogramSummary,
                     RunStatisticBarcode, RunSummaryBarcode, GroupRun, FlowcellHistogramSummary,
                     FlowcellStatisticBarcode,
                     FlowcellSummaryBarcode, Flowcell, Minion,
                     FlowcellTab, JobMaster, JobType)


class RunStatisticsBarcodeAdmin(admin.ModelAdmin):
    list_display = ('run_id', 'type')
    list_filter = ('run_id', )
    ordering = ('run_id', 'type__name')


class RunSummaryBarcodeAdmin(admin.ModelAdmin):
    list_display = ('run_id', 'barcode', 'type')
    list_filter = ('run_id', )
    ordering = ('run_id', 'barcode', 'type')


class FlowcellAdmin(admin.ModelAdmin):
    list_display = ('id', 'name', 'sample_name')


admin.site.register(FastqRead)
admin.site.register(FastqFile)
admin.site.register(FastqReadType)
admin.site.register(FlowcellStatisticBarcode)
admin.site.register(FlowcellSummaryBarcode)
admin.site.register(FlowcellHistogramSummary)
admin.site.register(FlowcellTab)
admin.site.register(Run)
admin.site.register(MinionEventType)
admin.site.register(MinionEvent)
admin.site.register(MinionScripts)
admin.site.register(MinionInfo)
admin.site.register(MinionRunInfo)
admin.site.register(MinionRunStats)
admin.site.register(MinionMessage)
admin.site.register(MinionControl)
admin.site.register(RunStatisticBarcode, RunStatisticsBarcodeAdmin)
admin.site.register(RunSummaryBarcode, RunSummaryBarcodeAdmin)
admin.site.register(ChannelSummary)
admin.site.register(HistogramSummary)
admin.site.register(UserOptions)
admin.site.register(Barcode)
admin.site.register(GroupRun)
admin.site.register(RunSummary)
admin.site.register(Flowcell, FlowcellAdmin)
admin.site.register(Minion)
admin.site.register(JobMaster)
admin.site.register(JobType)
