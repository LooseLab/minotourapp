# -*- coding: utf-8 -*-
# Generated by Django 1.11.6 on 2017-11-14 15:24
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('assembly', '0002_auto_20171113_1605'),
    ]

    operations = [
        migrations.AddField(
            model_name='gfasummary',
            name='allcontigs',
            field=models.TextField(default='[]'),
            preserve_default=False,
        ),
    ]
