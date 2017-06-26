# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-06-23 13:34
from __future__ import unicode_literals

from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0028_auto_20170623_1332'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='minioncontrol',
            name='setdate',
        ),
        migrations.AddField(
            model_name='minioncontrol',
            name='created_date',
            field=models.DateTimeField(auto_now_add=True, default=django.utils.timezone.now),
            preserve_default=False,
        ),
    ]
