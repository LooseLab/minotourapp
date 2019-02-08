from django.conf.urls import url

from communication import views

urlpatterns = [
    url(
        r'^api/v1/messages/$',
        views.new_messages_list,
        name="messages-list"),
    url(
        r'^messages/(?P<pk>[0-9]+)/$',
        views.message_details,
        name="message-detail"
    ),
]
