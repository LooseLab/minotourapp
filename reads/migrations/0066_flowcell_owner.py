# -*- coding: utf-8 -*-
# Generated by Django 1.11.6 on 2017-11-01 16:00
from __future__ import unicode_literals

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('reads', '0065_auto_20171101_1547'),
    ]

    operations = [
        migrations.AddField(
            model_name='flowcell',
            name='owner',
            field=models.ForeignKey(default=1, on_delete=django.db.models.deletion.CASCADE, related_name='flowcells', to=settings.AUTH_USER_MODEL),
            preserve_default=False,
        ),
    ]
