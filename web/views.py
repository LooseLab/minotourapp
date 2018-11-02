from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.forms import formset_factory
from django.http import JsonResponse
from django.shortcuts import redirect, render
from django.views.generic import ListView, DeleteView
from rest_framework.authtoken.models import Token

from communication.models import Message
from reads.models import Run, UserOptions, FastqRead, Experiment, Flowcell
from web.forms import SignUpForm, UserOptionsForm, ExperimentForm, ExperimentFlowcellForm

from django.contrib import messages


def index(request):

    if request.user.is_authenticated:

        return redirect('flowcells')

    else:

        return redirect('login')

def signup(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            raw_password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=raw_password)
            login(request, user)
            return redirect('flowcells')

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
def flowcells(request):
    flowcells = Flowcell.objects.filter(owner=request.user)
    return render(request, 'web/flowcells.html', context={'flowcells': flowcells})


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

    flowcell_id = int(request.GET.get('flowcell_id', 0))

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

        if search_value[0] == '>':

            reads_temp = FastqRead.objects\
                .filter(run__flowcell_id=flowcell_id)\
                .filter(run__flowcell__owner=request.user)\
                .filter(sequence_length__gt=search_value[1:])

        elif search_value[0] == '<':

            reads_temp = FastqRead.objects\
                .filter(run__flowcell_id=flowcell_id) \
                .filter(run__flowcell__owner=request.user)\
                .filter(sequence_length__lt=search_value[1:])

        else:

            reads_temp = FastqRead.objects \
                .filter(run__flowcell_id=flowcell_id) \
                .filter(run__flowcell__owner=request.user) \
                .filter(read_id__contains=search_value)

    else:

        reads_temp = FastqRead.objects\
            .filter(run__flowcell_id=flowcell_id)\
            .filter(run__flowcell__owner=request.user)\

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


class ExperimentList(ListView):
    model = Experiment


class ExperimentDelete(DeleteView):
    model = Experiment


@login_required
def experiments_create(request):

    ExperimentFlowcellFormSet = formset_factory(ExperimentFlowcellForm, extra=3)

    if request.method == 'POST':

        form = ExperimentForm(request.POST)
        experiment_flowcell_formset = ExperimentFlowcellFormSet(request.POST)

        if form.is_valid() and experiment_flowcell_formset.is_valid():

            experiment_name = form.cleaned_data['name']

            experiment = Experiment()
            experiment.name = experiment_name.upper()
            experiment.owner = request.user
            experiment.save()

            for experiment_flowcell_form in experiment_flowcell_formset:

                flowcell = experiment_flowcell_form.cleaned_data.get('flowcell')

                if flowcell:

                    flowcell.experiments.add(experiment)

            messages.success(request, 'The experiment {} was created with success.'.format(experiment.name))
            return redirect('experiment-list')

    else:

        form = ExperimentForm()
        experiment_flowcell_formset = ExperimentFlowcellFormSet()

    experiment_list = Experiment.objects.filter(owner=request.user)

    return render(
        request,
        'reads/experiments_create.html',
        {
            'form': form,
            'experiment_flowcell_formset': experiment_flowcell_formset,
            'experiment_list': experiment_list,
        }
    )

@login_required
def experiments_update(request, pk):

    experiment = Experiment.objects.get(pk=pk)

    ExperimentFlowcellFormSet = formset_factory(ExperimentFlowcellForm, extra=3)

    if request.method == 'POST':

        form = ExperimentForm(request.POST)
        experiment_flowcell_formset = ExperimentFlowcellFormSet(request.POST)

        if form.is_valid() and experiment_flowcell_formset.is_valid():

            experiment_name = form.cleaned_data['name']


            experiment.name = experiment_name.upper()
            experiment.save()

            for experiment_flowcell_form in experiment_flowcell_formset:

                flowcell = experiment_flowcell_form.cleaned_data.get('flowcell')

                if flowcell:

                    flowcell.experiments.add(experiment)

            messages.success(request, 'The experiment {} was updated with success.'.format(experiment.name))
            return redirect('experiment-list')

    else:

        form = ExperimentForm(instance=experiment)

        experiment_flowcell_list = experiment.flowcell_set.all()

        initial_data = []

        for experiment_flowcell in experiment_flowcell_list:

            initial_data.append({'flowcell': experiment_flowcell})

        experiment_flowcell_formset = ExperimentFlowcellFormSet(initial=initial_data)

    return render(
        request,
        'reads/experiments_update.html',
        {
            'form': form,
            'experiment_flowcell_formset': experiment_flowcell_formset,
        }
    )
