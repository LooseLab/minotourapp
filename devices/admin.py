from django.contrib import admin

from devices.models import Flowcell, MinION

admin.site.register(Flowcell)
admin.site.register(MinION)
