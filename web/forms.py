from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User


class UserOptionsForm(forms.Form):

    email = forms.EmailField(
        label='Email',
        max_length=120,
    )

    twitter_handle = forms.CharField(
        label='Twitter Handle',
        max_length=32,
        required=False,
    )

    receive_emails = forms.BooleanField(
        required=False,
    )

    receive_tweets = forms.BooleanField(
        required=False,
    )


class SignUpForm(UserCreationForm):

    first_name = forms.CharField(
        max_length=30,
        required=False,
        help_text="Optional."
    )

    last_name = forms.CharField(
        max_length=30,
        required=False,
        help_text='Optional.'
    )

    email = forms.EmailField(
        max_length=254,
        help_text='Required. Inform a valid email address.'
    )

    class Meta:
        model = User
        fields = (
            'username',
            'first_name',
            'last_name',
            'email',
            'password1',
            'password2',
        )
