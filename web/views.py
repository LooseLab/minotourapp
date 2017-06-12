from django.contrib.auth.decorators import login_required
from django.shortcuts import render
from reads.models import MinIONRun


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
    minion_runs = MinIONRun.objects.filter(owner=request.user)
    return render(request, 'web/previous_run.html', context={'minion_runs': minion_runs})

@login_required
def run_index(request, run_name):
    minion_runs = MinIONRun.objects.filter(run_name=run_name)
    return render(request, 'web/run_index.html', context={'minion_runs': minion_runs})
