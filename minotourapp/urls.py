from django.conf.urls import include
from django.conf.urls import url
from django.contrib import admin
#import notifications.urls
from rest_framework.authtoken import views as tok_views
from django.contrib.auth import views as auth_views

urlpatterns = [
    url(r'^login/$', auth_views.LoginView.as_view(template_name="registration/login.html"), name="login"),
    url(r'^logout/$', auth_views.LogoutView.as_view(), name="logout"),
    url(r'^change-password/$', auth_views.PasswordChangeView.as_view(template_name="registration/password_change.html"), name="password_change"),
    url(r'^change-password-done/$', auth_views.PasswordChangeDoneView.as_view(template_name="registration/password_change_done.html"), name="password_change_done"),
    url(r'^reset-password/$', auth_views.PasswordResetView.as_view(template_name="registration/password_reset.html"), name="password_reset"),
    url(r'^reset-password-done/$', auth_views.PasswordResetDoneView.as_view(template_name="registration/password_reset_done.html"), name="password_reset_done"),
    url(r'^reset/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$', auth_views.PasswordResetConfirmView.as_view(template_name="registration/password_reset_confirm.html"), name='password_reset_confirm'),
    url(r'^admin/', admin.site.urls),
    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    url(r'^web/', include('web.urls')),
    url(r'^', include('reads.urls')),
    url(r'^', include('alignment.urls')),
    url(r'^', include('reference.urls')),
    url(r'^', include('communication.urls')),
    url(r'^', include('minikraken.urls')),
    url(r'^$', auth_views.LoginView.as_view(template_name="registration/login.html")),
    url(r'^api-token-auth/', tok_views.obtain_auth_token),
#    url('^inbox/notifications/', include(notifications.urls, namespace='notifications')),
]
