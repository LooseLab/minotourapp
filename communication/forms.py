"""
Defines form handling for creating notification conditio
"""
from django import forms

from .models import NotificationConditions


class NotificationConditionForm(forms.ModelForm):
    """
    Form to save references uploaded by the user
    """
    class Meta:
        model = NotificationConditions
        fields = ("upper_limit", "lower_limit", "notification_type", "flowcell", "chromosome", "reference", "coverage_target", "run_until")
