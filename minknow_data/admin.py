from django.contrib import admin

from minknow_data.models import Minion, MinionControl, MinionInfo, MinionRunStats, MinionRunInfo, MinionEventType, \
    MinionEvent, MinionScripts, MinionMessage

# Register your models here.

admin.site.register(MinionEventType)
admin.site.register(MinionEvent)
admin.site.register(MinionScripts)
admin.site.register(MinionInfo)
admin.site.register(MinionRunInfo)
admin.site.register(MinionRunStats)
admin.site.register(MinionMessage)
admin.site.register(MinionControl)
admin.site.register(Minion)
