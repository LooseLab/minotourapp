# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-06-26 14:09
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0031_minionrunstatus_minknow_version'),
    ]

    operations = [
        migrations.AddField(
            model_name='minionstatus',
            name='minKNOW_disk_available',
            field=models.FloatField(blank=True, null=True),
        ),
    ]
