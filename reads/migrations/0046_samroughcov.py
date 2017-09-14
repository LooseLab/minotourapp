# -*- coding: utf-8 -*-
# Generated by Django 1.11.1 on 2017-09-13 15:45
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0045_runstatisticbarcode_channel_presence'),
    ]

    operations = [
        migrations.CreateModel(
            name='SamRoughCov',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reference', models.TextField()),
                ('chromosome', models.TextField()),
                ('position', models.IntegerField()),
                ('incdel', models.IntegerField()),
                ('read_type', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='rc_type', to='reads.FastqReadType')),
                ('run_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='rc_runid', to='reads.MinIONRun')),
            ],
        ),
    ]
