# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-06-16 19:12
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0020_minionstatus_run_id'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='minionstatus',
            name='run_id',
        ),
        migrations.AddField(
            model_name='minionrunstats',
            name='minKNOW_read_count',
            field=models.IntegerField(default=0),
        ),
    ]
