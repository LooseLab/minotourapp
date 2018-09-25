from django.contrib import admin

from reads.models import Flowcell, MinION

admin.site.register(Flowcell)
admin.site.register(MinION)
