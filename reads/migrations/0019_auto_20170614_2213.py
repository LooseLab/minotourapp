# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-06-14 22:13
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0018_auto_20170614_2059'),
    ]

    operations = [
        migrations.AlterField(
            model_name='minionrun',
            name='minION',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='minionrun', to='reads.MinION'),
        ),
    ]
