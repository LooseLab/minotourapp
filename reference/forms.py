from django import forms

from .models import ReferenceInfo


class ReferenceForm(forms.ModelForm):
    """
    Form to save references uploaded by the user
    """
    class Meta:
        model = ReferenceInfo
        fields = ('file_location', 'file_name')
