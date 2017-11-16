from django.contrib import admin

from .models import GfaStore
from .models import GfaSummary

#Register your models here.

admin.site.register(GfaStore)
admin.site.register(GfaSummary)
