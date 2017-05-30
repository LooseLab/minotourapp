from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from reads.models import MinionRun


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
def external_links(request):
    return render(request, 'web/external_links.html')


@login_required
def minup(request):
    return render(request, 'web/minup.html')


@login_required
def tutorial(request):
    return render(request, 'web/tutorial.html')


@login_required
def previous_run(request):
    minion_runs = MinionRun.objects.filter(owner=request.user)
    return render(request, 'web/previous_run.html', context={'minion_runs': minion_runs})
