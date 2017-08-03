from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from reads.models import MinIONRun
from django.db.models import Q
from datetime import datetime, timedelta
from django.utils import timezone

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