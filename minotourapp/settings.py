"""
Django settings for minotourapp project.

Generated by 'django-admin startproject' using Django 1.11.1.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.11/ref/settings/
"""
import os
import tempfile

from celery.schedules import crontab
from kombu import Exchange, Queue

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
from minotourapp.utils import get_env_variable

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = get_env_variable("MT_SECRET_KEY")

# DEBUG = bool(os.environ.get('MT_DJANGO_DEBUG', True))
ALLOWED_HOSTS = ['*', ]

# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.humanize',
    'rest_framework',
    'rest_framework.authtoken',
    'django_extensions',
    'reads',
    'web',
    'alignment',
    'reference',
    'communication',
    'metagenomics',
    'readuntil',
    'guardian',
    'artic',
    'minknow_data'
]

MIDDLEWARE = [
    'minotourapp.utils.GZIPBodyDecompressMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'whitenoise.middleware.WhiteNoiseMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'minotourapp.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, 'web/templates/')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'minotourapp.wsgi.application'

CACHES = {
    "default": {
        "BACKEND": "django_redis.cache.RedisCache",
        "LOCATION": "redis://127.0.0.1:6379/1",
        "OPTIONS": {
            "CLIENT_CLASS": "django_redis.client.DefaultClient",
        }
    }
}

DATABASES = {
    'default': {
        'ENGINE': get_env_variable("MT_DB_ENGINE"),
        'USER': get_env_variable("MT_DB_USER"),
        'PASSWORD': get_env_variable("MT_DB_PASS"),
        'HOST': get_env_variable("MT_DB_HOST"),
        'PORT': get_env_variable("MT_DB_PORT"),
        'NAME': get_env_variable("MT_DB_NAME"),
        'CONN_MAX_AGE': 0,
    }
}

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/

LANGUAGE_CODE = 'en-gb'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)

STATIC_URL = '/static/'
STATIC_ROOT = os.path.join(BASE_DIR, '../minotour_static/')
MEDIA_ROOT = os.path.join(BASE_DIR, '../data/references')
MEDIA_URL = '/media/'

REST_FRAMEWORK = {
    'DEFAULT_PERMISSION_CLASSES': [
        # 'rest_framework.permissions.DjangoModelPermissionsOrAnonReadOnly'
        'rest_framework.permissions.IsAuthenticated'
    ],
    'DEFAULT_AUTHENTICATION_CLASSES': (
        'rest_framework.authentication.BasicAuthentication',
        'rest_framework.authentication.SessionAuthentication',
        'rest_framework.authentication.TokenAuthentication',
    ),
    # 'DEFAULT_RENDERER_CLASSES': [
    #     'minotourapp.renderers.GzipJsonRenderer',
    #     'rest_framework.renderers.BrowsableAPIRenderer',
    # ]
}

LOGIN_URL = '/login'
LOGIN_REDIRECT_URL = '/web/private/flowcells'
LOGOUT_REDIRECT_URL = LOGIN_URL

# EMAIL_BACKEND = 'django_mailgun.MailgunBackend'
EMAIL_BACKEND = 'django.core.mail.backends.console.EmailBackend'
MAILGUN_ACCESS_KEY = get_env_variable("MT_MAILGUN_ACCESS_KEY")
MAILGUN_SERVER_NAME = get_env_variable("MT_MAILGUN_SERVER_NAME")

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '{levelname} {asctime} {module} {process:d} {thread:d} {message}',
            'style': '{',
        },
    },
    'handlers': {
        'file': {
            'level': 'ERROR',
            'class': 'logging.FileHandler',
            'filename': os.path.join(BASE_DIR, 'minotour.log'),
        },
        'console': {
            'level': 'INFO',
            'class': 'logging.StreamHandler',
            'formatter': 'verbose',
        }
    },
    'loggers': {
        'django': {
            'handlers': ['file', 'console'],
            'level': 'INFO',
            'propagate': True,
        },
    },
}

default_exchange = Exchange('default', type='direct')
priority_exchange = Exchange('priority_queue', type='direct')
CELERY_QUEUES = (
    Queue('minimap', default_exchange, routing_key='alignment#', consumer_arguments={'x-priority': 5}),
)
CELERY_IMPORTS = ('web.tasks', 'communication.tasks_send_message',
                  'alignment.tasks_alignment', 'reads.tasks.tasks_archive_flowcell')
CELERY_ROUTES = ({
    'run_minimap2_alignment': {
        'queue': 'minimap',
        'routing_key': "alignment#"
    },
    'remove_old_references': {
        'queue': 'minimap',
        'routing_key': "alignment#"
    },
    'web.tasks_centrifuge': {
        'queue': 'metagenomics',
        'routing_key': 'metagenomics'
    }
})

CELERY_BROKER_URL = get_env_variable("MT_CELERY_BROKER_URL")
CELERY_RESULT_BACKEND = get_env_variable("MT_CELERY_RESULT_BACKEND")
# Celery Data Format
CELERY_ACCEPT_CONTENT = ['application/json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_RESULT_SERIALIZER = 'json'
CELERY_TIMEZONE = 'Europe/London'
CELERY_BEAT_SCHEDULE = {
    'rapid-monitor': {
        'task': 'web.tasks.run_monitor',
        'schedule': 30,
    },
    'send_messages': {
        'task': 'communication.tasks_send_message.send_messages',
        'schedule': 30,
    },
    'check_condition_is_met': {
        'task': 'communication.tasks_send_message.check_condition_is_met',
        'schedule': 30,
    },
    'remove_old_references': {
        'task': 'alignment.tasks_alignment.remove_old_references',
        'schedule': 1800,
        'options': {
            "queue": "minimap"
        }
    },
    # Run this every midnight
    'secure_artic_runs': {
        'task': 'artic.task_secure_artic_runs.secure_artic_runs',
        'schedule': crontab(minute=0, hour=0)
    },
    # Run the thing on the thing
    'clear_intermediate_alignment': {
        'task': 'alignment.tasks_alignment.aggregate_intermediate_table',
        'schedule': 45,
        'options': {
            "queue": "minimap"
        }
    }
}

# For sending twitter messages
TWITTOKEN = get_env_variable("MT_TWITTOKEN")
TWITTOKEN_SECRET = get_env_variable("MT_TWITTOKEN_SECRET")
TWITCONSUMER_KEY = get_env_variable("MT_TWITCONSUMER_KEY")
TWITCONSUMER_SECRET = get_env_variable("MT_TWITCONSUMER_SECRET")

# Variables for storing additonal files
# References
# REFERENCE_LOCATION = get_env_variable("MT_REFERENCE_LOCATION")
#
# if REFERENCE_LOCATION == '':
#     REFERENCE_LOCATION = os.path.join(BASE_DIR, 'data')

PAGINATION_PAGE_SIZE = 1000

MINIMAP2 = get_env_variable("MT_MINIMAP2")

USE_X_FORWARDED_HOST = True

# SECURITY WARNING: don't run with debug turned on in production!

DEBUG = bool(int(get_env_variable("MT_DJANGO_DEBUG")))

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',  # default
    'guardian.backends.ObjectPermissionBackend',
)

# if DEBUG:
#    INTERNAL_IPS = ('127.0.0.1', 'localhost',)
#    MIDDLEWARE += (
#        'debug_toolbar.middleware.DebugToolbarMiddleware',
#    )
#
#    INSTALLED_APPS += (
#        'debug_toolbar',
#    )
#
#    DEBUG_TOOLBAR_PANELS = [
#        'debug_toolbar.panels.versions.VersionsPanel',
#        'debug_toolbar.panels.timer.TimerPanel',
#        'debug_toolbar.panels.settings.SettingsPanel',
#        'debug_toolbar.panels.headers.HeadersPanel',
#        'debug_toolbar.panels.request.RequestPanel',
#        'debug_toolbar.panels.sql.SQLPanel',
#        'debug_toolbar.panels.staticfiles.StaticFilesPanel',
#        'debug_toolbar.panels.templates.TemplatesPanel',
#        'debug_toolbar.panels.cache.CachePanel',
#        'debug_toolbar.panels.signals.SignalsPanel',
#        'debug_toolbar.panels.logging.LoggingPanel',
#        'debug_toolbar.panels.redirects.RedirectsPanel',
#    ]
#
#    DEBUG_TOOLBAR_CONFIG = {
#        'INTERCEPT_REDIRECTS': False,
#    }
