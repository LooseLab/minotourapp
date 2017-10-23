# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.contrib.auth.hashers import make_password
from django.contrib.auth.models import User
from django.db import migrations

from rest_framework.authtoken.models import Token

def forwards_func(apps, schema_editor):
    # We get the model from the versioned app registry;
    # if we directly import it, it'll be the wrong version
    db_alias = schema_editor.connection.alias

    FastqReadType = apps.get_model("reads", "FastqReadType")

    FastqReadType.objects.using(db_alias).bulk_create([
        FastqReadType(name="Template"),
        FastqReadType(name="Complement"),
        FastqReadType(name="2D"),
        FastqReadType(name="1D^2")
    ])

    JobType = apps.get_model("reads", "JobType")

    JobType.objects.using(db_alias).bulk_create([
        JobType(name="ChanCalc", private=True),
        JobType(name="ProcAlign", private=True),
        JobType(name="Alignment", private=True),
        JobType(name="Minimap2", private=True)
    ])

    MinIONEventType = apps.get_model("reads", "MinIONEventType")

    MinIONEventType.objects.using(db_alias).bulk_create([
        MinIONEventType(name="active"),
        MinIONEventType(name="inactive"),
        MinIONEventType(name="sequencing"),
        MinIONEventType(name="unplugged"),
        MinIONEventType(name="connected"),
    ])

    User.objects.using(db_alias).bulk_create([
        User(
            username="demo",
            is_active=True,
            is_superuser=True,
            is_staff=True,
            email="demo@minotour.org"
        )
    ])

    users = User.objects.filter(username="demo")

    for user in users:
        user.password = make_password('asdf1234')
        user.save()

        Token.objects.get_or_create(user=user)


def reverse_func(apps, schema_editor):
    # forwards_func() creates two Country instances,
    # so reverse_func() should delete them.

    db_alias = schema_editor.connection.alias

    FastqReadType = apps.get_model("reads", "FastqReadType")

    FastqReadType.objects.using(db_alias).filter(name="Template").delete()
    FastqReadType.objects.using(db_alias).filter(name="Complement").delete()
    FastqReadType.objects.using(db_alias).filter(name="2D").delete()
    FastqReadType.objects.using(db_alias).filter(name="1D^2").delete()

    JobType = apps.get_model("reads", "JobType")

    JobType.objects.using(db_alias).filter(name="ChanCalc").delete()
    JobType.objects.using(db_alias).filter(name="ProcAlign").delete()
    JobType.objects.using(db_alias).filter(name="Alignment").delete()
    JobType.objects.using(db_alias).filter(name="Minimap2").delete()

    MinIONEventType = apps.get_model("reads", "MinIONEventType")

    MinIONEventType.objects.using(db_alias).filter(name="active").delete()
    MinIONEventType.objects.using(db_alias).filter(name="inactive").delete()
    MinIONEventType.objects.using(db_alias).filter(name="sequencing").delete()
    MinIONEventType.objects.using(db_alias).filter(name="unplugged").delete()
    MinIONEventType.objects.using(db_alias).filter(name="connected").delete()

    User.objects.using(db_alias).filter(username="demo").delete()


class Migration(migrations.Migration):

    dependencies = [
        ('reads', '0054_minionrun_to_delete'),
    ]

    operations = [
        migrations.RunPython(forwards_func, reverse_func),
    ]
