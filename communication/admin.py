from django.contrib import admin

from communication.models import Message,NotificationConditions

admin.site.register(Message)
admin.site.register(NotificationConditions)
