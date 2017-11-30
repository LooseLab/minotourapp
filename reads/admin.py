from django.contrib import admin

from .models import FastqRead
from .models import RunSummaryBarcode
from .models import RunStatisticBarcode
from .models import FastqReadType
from .models import MinION
from .models import MinIONControl
from .models import MinIONEvent
from .models import MinIONEventType
from .models import MinIONRun
from .models import MinIONRunStats
from .models import MinIONRunStatus
from .models import MinIONScripts
from .models import MinIONStatus
from .models import MinIONmessages
from .models import JobMaster
from .models import JobType
from .models import ChannelSummary
from .models import HistogramSummary
from .models import UserOptions
from .models import Barcode
from .models import BarcodeGroup
from .models import FlowCell
from .models import FlowCellRun
from .models import GroupRunType


class MinionRunAdmin(admin.ModelAdmin):
    list_display = ('start_time', 'active','sample_name', 'minKNOW_version', 'minKNOW_flow_cell_id', 'run_name', 'run_id', 'is_barcoded','minION', 'id', 'last_read', 'last_entry')


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
admin.site.register(MinIONRun, MinionRunAdmin)
admin.site.register(MinION)
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
admin.site.register(FlowCell)
admin.site.register(FlowCellRun)
admin.site.register(GroupRunType)
