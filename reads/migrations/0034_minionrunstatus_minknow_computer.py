# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-08-28 19:29
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        #('reads', '0033_auto_20170808_1446'),
        ('reads', '0034_auto_20170824_1127'),
    ]

    operations = [
        migrations.AddField(
            model_name='minionrunstatus',
            name='minKNOW_computer',
            field=models.TextField(blank=True, max_length=128, null=True),
        ),
    ]
