from django.contrib import admin

from .models import FastqRead
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
admin.site.register(RunSummary)
