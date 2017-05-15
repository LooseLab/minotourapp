from django.contrib import admin

from .models import FastqRead
from .models import FastqReadType
from .models import MinionRun


admin.site.register(FastqRead)
admin.site.register(FastqReadType)
admin.site.register(MinionRun)
