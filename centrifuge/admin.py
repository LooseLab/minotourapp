from django.contrib import admin

from centrifuge.models import CentOutput, LineageKey, LineageValues

admin.site.register(CentOutput)
admin.site.register(LineageKey)
admin.site.register(LineageValues)
