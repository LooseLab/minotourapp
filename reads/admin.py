from django.contrib import admin

from .models import FastqRead, RunSummaryBarCode
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
from .models import RunStatistic
from .models import RunSummary


class RunSummaryAdmin(admin.ModelAdmin):
    list_display = ('run_id', 'type')
    list_filter = ('run_id', )
    ordering = ('run_id', 'type__name')

class RunSummaryBarcodeAdmin(admin.ModelAdmin):
    list_display = ('run_id', 'barcode', 'type')
    list_filter = ('run_id', )
    ordering = ('run_id', 'barcode', 'type')

admin.site.register(FastqRead)
admin.site.register(FastqReadType)
admin.site.register(MinIONRun)
admin.site.register(RunStatistic)
admin.site.register(MinION)
admin.site.register(MinIONEventType)
admin.site.register(MinIONEvent)
admin.site.register(MinIONScripts)
admin.site.register(MinIONStatus)
admin.site.register(MinIONRunStatus)
admin.site.register(MinIONRunStats)
admin.site.register(MinIONmessages)
admin.site.register(MinIONControl)
admin.site.register(RunSummary, RunSummaryAdmin)
admin.site.register(RunSummaryBarCode, RunSummaryBarcodeAdmin)

