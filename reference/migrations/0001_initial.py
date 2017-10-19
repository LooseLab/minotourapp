# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-14 15:06
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='ReferenceInfo',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reference_name', models.CharField(max_length=512)),
                ('filename', models.CharField(max_length=512)),
                ('bwa_index_file_location', models.CharField(max_length=512, null=True)),
                ('minimap2_index_file_location', models.CharField(max_length=512, null=True)),
                ('totalrefleN', models.IntegerField()),
                ('private', models.BooleanField(default=False)),
            ],
        ),
        migrations.CreateModel(
            name='ReferenceLines',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('line_name', models.CharField(max_length=512)),
                ('chromosome_length', models.IntegerField()),
                ('reference', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='referencelines', to='reference.ReferenceInfo')),
            ],
        ),
    ]
