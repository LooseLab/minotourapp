# -*- coding: utf-8 -*-
# Generated by Django 1.11.4 on 2017-10-10 09:56
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0052_auto_20171010_0941'),
    ]

    operations = [
        migrations.RenameField(
            model_name='jobmaster',
            old_name='job_name',
            new_name='job_type',
        ),
        migrations.RenameField(
            model_name='jobmaster',
            old_name='run_id',
            new_name='run',
        ),
    ]
