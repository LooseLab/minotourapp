# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-13 17:32
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0046_samroughcov'),
    ]

    operations = [
        migrations.CreateModel(
            name='PafStore',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('qsn', models.CharField(max_length=256)),
                ('qsl', models.IntegerField()),
                ('qs', models.IntegerField()),
                ('qe', models.IntegerField()),
                ('rs', models.CharField(max_length=1)),
                ('tsn', models.CharField(max_length=256)),
                ('tsl', models.IntegerField()),
                ('ts', models.IntegerField()),
                ('te', models.IntegerField()),
                ('nrm', models.IntegerField()),
                ('abl', models.IntegerField()),
                ('mq', models.IntegerField()),
                ('read_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='pafreadalignment', to='reads.FastqRead')),
                ('run_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='pafalignemnts', to='reads.MinIONRun')),
            ],
        ),
    ]
