from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.http import JsonResponse
from django.shortcuts import redirect, render
from rest_framework.authtoken.models import Token

from communication.models import Message
from devices.models import Flowcell
from reads.models import Run, UserOptions, FastqRead
from web.forms import SignUpForm, UserOptionsForm


def signup(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            raw_password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=raw_password)
            login(request, user)
            return redirect('private-index')

    else:
        form = SignUpForm()

    return render(request, 'web/signup.html', {'form': form})


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
def minup(request):
    return render(request, 'web/minup.html')


@login_required
def tutorial(request):
    return render(request, 'web/tutorial.html')


@login_required
def runs(request):
    minion_runs = Run.objects.filter(owner=request.user)
    return render(request, 'web/runs.html', context={'minion_runs': minion_runs})


@login_required
def flowcells(request):
    flowcells = Flowcell.objects.filter(owner=request.user)
    return render(request, 'web/flowcells.html', context={'flowcells': flowcells})


@login_required
def run_index(request, pk):
    minion_run = Run.objects.get(pk=pk)
    return render(request, 'web/run_index.html', context={'minion_run': minion_run})


@login_required
def flowcell_index(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)
    return render(request, 'web/flowcell_index.html', context={'flowcell': flowcell})


@login_required
def flowcell_reads(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    return render(request, 'web/flowcell_reads.html', context={'flowcell': flowcell})


def flowcell_reads_data(request):

    # column_0_data = request.GET.get('column[0][data]', '')
    # column_0_name = request.GET.get('column[0][name]', '')
    # column_0_searchable = request.GET.get('column[0][searchable]', '')
    # column_0_orderable = request.GET.get('column[0][orderable]', '')
    # column_0_search_value = request.GET.get('column[0][search][value]', '')
    # column_0_search_regex = request.GET.get('column[0][search][regex]', '')

    query_columns = [
        'read_id',
        'read',
        'channel',
        'sequence_length',
        'run__runid',
        'barcode__name',
    ]

    query_columns_string = ['read_id', 'read', 'channel', 'sequence_length', 'run__runid', 'barcode__name']

    draw = int(request.GET.get('draw', 0))

    search_value = request.GET.get('search[value]', '')

    start = int(request.GET.get('start', 0))

    length = int(request.GET.get('length', 10))

    end = start + length
    # Which column s
    order_column = request.GET.get('order[0][column]', '')
    # ascending descending
    order_dir = request.GET.get('order[0][dir]', '')

    if not search_value == "":
        reads_temp = FastqRead.objects.filter(read_id__contains=search_value)

    else:
        reads_temp = FastqRead.objects.all()

    if order_column:

        if order_dir == 'desc':

            reads_temp2 = reads_temp.order_by('-{}'.format(query_columns[int(order_column)]))

        else:

            reads_temp2 = reads_temp.order_by('{}'.format(query_columns[int(order_column)]))

    reads = reads_temp2.values('run__runid', 'barcode__name', 'read_id', 'read', 'channel', 'sequence_length')

    records_total = len(reads)

    result = {
        'draw': draw,
        "recordsTotal": records_total,
        "recordsFiltered": records_total,
        "data": list(reads[start:end])
    }

    return JsonResponse(result, safe=True)


@login_required
def remotecontrol(request):
    return render(request, 'web/remotecontrol.html')


@login_required
def sandbox(request):
    return render(request, 'web/sandbox.html')
