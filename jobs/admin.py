from django.contrib import admin

from jobs.models import JobType, JobMaster

admin.site.register(JobMaster)
admin.site.register(JobType)