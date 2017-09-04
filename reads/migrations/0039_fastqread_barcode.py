# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-04 11:29
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0038_remove_fastqread_barcode'),
    ]

    operations = [
        migrations.AddField(
            model_name='fastqread',
            name='barcode',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='reads', to='reads.Barcode'),
        ),
    ]
