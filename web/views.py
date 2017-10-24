from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.shortcuts import render
from rest_framework.authtoken.models import Token

from communication.models import Message
from reads.models import MinIONRun
from reads.models import UserOptions
from web.forms import UserOptionsForm


def index(request):
    return render(request, 'web/index.html')


def current(request):
    return render(request, 'web/current_run.html')


def log_in(request):
    return render(request, 'web/log_in.html')


@login_required
def private_index(request):
    return render(request, 'web/private_index.html')


@login_required
def profile(request):

    try:

        user_options = UserOptions.objects.get(owner=request.user)

    except ObjectDoesNotExist:

        user_options = UserOptions.objects.create(owner=request.user)

    auth_token = Token.objects.get(user=request.user)

    messages = Message.objects.filter(recipient=request.user).order_by('-created_date')

    if request.method == 'POST':

        form = UserOptionsForm(request.POST)

        if form.is_valid():

            user = User.objects.get(pk=user_options.owner_id)
            user.email = form.cleaned_data['email']
            user.save()

            user_options.twitterhandle = form.cleaned_data['twitter_handle']
            user_options.tweet = form.cleaned_data['receive_tweets']
            user_options.email = form.cleaned_data['receive_emails']
            user_options.save()

    else:

        form = UserOptionsForm()

    return render(
        request, 'web/profile.html',
        context={
            'authToken': auth_token,
            'userDetails': user_options,
            'messages': messages,
            'form': form,
        }
    )


@login_required
def external_links(request):
    return render(request, 'web/external_links.html')


@login_required
def minup(request):
    return render(request, 'web/minup.html')


@login_required
def tutorial(request):
    return render(request, 'web/tutorial.html')


@login_required
def current_run(request):
    return render(request, 'web/current_run2.html')


@login_required
def previous_run(request):
    minion_runs = MinIONRun.objects.filter(owner=request.user)
    return render(request, 'web/previous_run.html', context={'minion_runs': minion_runs})


@login_required
def run_index(request, pk):
    minion_run = MinIONRun.objects.get(pk=pk)
    return render(request, 'web/run_index.html', context={'minion_run': minion_run})


@login_required
def remotecontrol(request):
    return render(request, 'web/remotecontrol.html')


@login_required
def prevremotecontrol(request):
    return render(request, 'web/prevremotecontrol.html')
