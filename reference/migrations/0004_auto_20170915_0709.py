# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-15 07:09
from __future__ import unicode_literals

from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('reference', '0003_auto_20170914_2113'),
    ]

    operations = [
        migrations.AddField(
            model_name='referenceinfo',
            name='users',
            field=models.ManyToManyField(blank=True, to=settings.AUTH_USER_MODEL),
        ),
    ]
