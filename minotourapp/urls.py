from django.conf.urls import include, url
from django.contrib import admin
from django.contrib.auth import views as auth_views
from rest_framework.authtoken import views as tok_views


urlpatterns = [
    url(r'^login/$', auth_views.LoginView.as_view(template_name="registration/login.html"), name="login",
        kwargs={'redirect_authenticated_user': True}),

    url(r'^logout/$', auth_views.LogoutView.as_view(), name="logout"),

    url(r'^change-password/$', auth_views.PasswordChangeView.as_view(template_name="registration/password_change.html"), name="password_change"),

    url(r'^change-password-done/$', auth_views.PasswordChangeDoneView.as_view(template_name="registration/password_change_done.html"), name="password_change_done"),

    url(r'^password-reset/$', auth_views.password_reset, name="password_reset"),

    url(r'^password-reset-done/$', auth_views.password_reset_done, name="password_reset_done"),

    url(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$', auth_views.password_reset_confirm, name='password_reset_confirm'),

    url(r'^reset/done/$', auth_views.password_reset_complete, name='password_reset_complete'),

    url(r'^admin/', admin.site.urls),

    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),

    url(r'^web/', include('web.urls')),
    url(r'^', include('reads.urls')),
    url(r'^', include('alignment.urls')),
    url(r'^', include('assembly.urls')),
    url(r'^', include('reference.urls')),
    url(r'^', include('communication.urls')),
    url(r'^', include('minikraken.urls')),
#    url(r'^', include('tabs.urls')),
    url(r'^$', auth_views.LoginView.as_view(template_name="registration/login.html")),
    url(r'^api-token-auth/', tok_views.obtain_auth_token),
#    url('^inbox/notifications/', include(notifications.urls, namespace='notifications')),
]
