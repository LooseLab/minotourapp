from django.contrib import admin

from .models import FastqRead
from .models import FastqReadType
from .models import MinIONRun
from .models import RunStatistic
from .models import MinION
from .models import MinIONEvent
from .models import MinIONEventType

admin.site.register(FastqRead)
admin.site.register(FastqReadType)
admin.site.register(MinIONRun)
admin.site.register(RunStatistic)
admin.site.register(MinION)
admin.site.register(MinIONEventType)
admin.site.register(MinIONEvent)