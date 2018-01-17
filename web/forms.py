from django import forms


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

