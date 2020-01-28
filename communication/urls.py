from django.conf.urls import url
from communication import views

urlpatterns = [
    url(
        r'^api/v1/messages/create_conditions$',
        views.get_or_create_conditions,
        name="create_conditions"
    ),
    url(
        r'^api/v1/messages/retrieve_messages',
        views.retrieve_messages,
        name="retrieve_messages"
    ),
    url(
        r'^api/v1/messages/retrieve_conditions$',
        views.get_or_create_conditions,
        name="retrieve_conditions"
    ),
]
