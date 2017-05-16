from django.contrib import admin

from .models import FastqRead
from .models import FastqReadType
from .models import MinionRun
from .models import RunStatistic

admin.site.register(FastqRead)
admin.site.register(FastqReadType)
admin.site.register(MinionRun)
admin.site.register(RunStatistic)
