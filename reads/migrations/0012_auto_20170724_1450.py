# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-07-24 14:50
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0011_sam'),
    ]

    operations = [
        migrations.AlterField(
            model_name='sam',
            name='cigar',
            field=models.TextField(),
        ),
        migrations.AlterField(
            model_name='sam',
            name='qual',
            field=models.TextField(),
        ),
        migrations.AlterField(
            model_name='sam',
            name='seq',
            field=models.TextField(),
        ),
    ]