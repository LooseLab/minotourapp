# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-10-01 20:31
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('minikraken', '0005_auto_20171001_1953'),
    ]

    operations = [
        migrations.AlterField(
            model_name='parsedkraken',
            name='NCBItaxid',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='parsedkraken',
            name='directreads',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='parsedkraken',
            name='indentation',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='parsedkraken',
            name='percentage',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='parsedkraken',
            name='rank',
            field=models.CharField(blank=True, max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='parsedkraken',
            name='rootreads',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='parsedkraken',
            name='sci_name',
            field=models.CharField(blank=True, max_length=128, null=True),
        ),
    ]
