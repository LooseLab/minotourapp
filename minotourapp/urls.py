from django.conf.urls import include, url
from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.urls import path
from rest_framework.authtoken import views as tok_views

from web.views import index

urlpatterns = [
    url(r'^login/$', auth_views.LoginView.as_view(template_name="registration/login.html"), name="login",
        kwargs={'redirect_authenticated_user': True}),

    url(r'^logout/$', auth_views.LogoutView.as_view(), name="logout"),

    url(r'^change-password/$', auth_views.PasswordChangeView.as_view(template_name="registration/password_change.html"), name="password_change"),

    url(r'^change-password-done/$', auth_views.PasswordChangeDoneView.as_view(template_name="registration/password_change_done.html"), name="password_change_done"),

    url(r'^admin/', admin.site.urls),

    url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),

    url(r'^web/', include('web.urls')),
    url(r'^', include('reads.urls')),
    url(r'^', include('alignment.urls')),
    url(r'^', include('reference.urls')),
    url(r'^', include('communication.urls')),
    url(r'^', include('metagenomics.urls')),
    url(r'^', include('artic.urls')),
    # url(r'^', include('readuntil.urls')),
    url(r'^', include('minknow_data.urls')),
    url(r'^$', index),
    url(r'^api-token-auth/', tok_views.obtain_auth_token),
    path('accounts/', include('django.contrib.auth.urls')),
]
