# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-04 11:41
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0040_auto_20170904_1138'),
    ]

    operations = [
        migrations.AddField(
            model_name='runstatisticbarcode',
            name='barcode',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='runstatistics', to='reads.Barcode'),
        ),
        migrations.AddField(
            model_name='runsummarybarcode',
            name='barcode',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='runsummaries', to='reads.Barcode'),
        ),
    ]