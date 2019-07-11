from django.conf.urls import url

from jobs import views

urlpatterns = [
    url(
        r'^api/v1/tasks/reset$',
        views.reset_task,
        name="restart-task",
    ),
    url(
        r'^api/v1/tasktypes/$',
        views.task_types_list,
        name="task-types-list"),
    # url(  # to be refactored
    #     r'^api/v1/runs/(?P<pk>[0-9]+)/settask/$',
    #     views.set_task_detail_all,
    #     name="set-task-detail-all"),
    url(
        r'^api/v1/tasks/$',
        views.get_or_create_tasks,
        name='task_list'
    )
]
