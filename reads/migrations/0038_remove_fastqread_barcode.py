# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-04 11:23
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0037_merge_20170904_1001'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='fastqread',
            name='barcode',
        ),
    ]
