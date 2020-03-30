from django.conf.urls import url
from communication import views

urlpatterns = [
    url(
        r'^api/v1/messages/conditions$',
        views.get_create_delete_conditions,
        name="conditions"
    ),
]
