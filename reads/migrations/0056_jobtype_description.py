# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-10-24 18:15
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0056_auto_20171026_1235'),
    ]

    operations = [
        migrations.AddField(
            model_name='jobtype',
            name='description',
            field=models.TextField(blank=True, null=True),
        ),
    ]
