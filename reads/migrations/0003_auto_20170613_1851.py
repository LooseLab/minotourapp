# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-06-13 18:51
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0002_auto_20170613_1423'),
    ]

    operations = [
        migrations.AddField(
            model_name='minionrunstats',
            name='asic_temp',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='minionrunstats',
            name='heat_sink_temp',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='minionrunstats',
            name='in_strand',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='minionrunstats',
            name='mean_ratio',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='minionrunstats',
            name='open_pore',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='minionrunstats',
            name='voltage_value',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='minionrunstats',
            name='event_yield',
            field=models.IntegerField(blank=True, null=True),
        ),
    ]
