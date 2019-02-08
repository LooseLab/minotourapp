from django.conf.urls import url

from jobs import views

urlpatterns = [
    url(
        r'^api/v1/tasktypes/$',
        views.task_types_list,
        name="task-types-list"),
    url(  # TODO DELETE
        r'^api/v1/runs/(?P<pk>[0-9]+)/tasks/$',
        views.tasks_detail_all,
        name="tasks-detail-all"),
    url(  # to be refactored
        r'^api/v1/runs/(?P<pk>[0-9]+)/settask/$',
        views.set_task_detail_all,
        name="set-task-detail-all"),
    url(
        r'^api/v1/tasks/$',
        views.task_list,
        name='task_list'
    )
]
