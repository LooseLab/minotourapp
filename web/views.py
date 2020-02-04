import pandas as pd
from django.contrib import messages
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.forms import formset_factory
from django.http import HttpResponse, JsonResponse
from django.shortcuts import redirect, render
from django.views.generic import ListView, DeleteView
from rest_framework.authtoken.models import Token

from centrifuge.models import CentrifugeOutput
from communication.models import Message
from reads.models import Run, UserOptions, FastqRead, Experiment, Flowcell, MinIONRunStats, JobType, JobMaster
from web.forms import SignUpForm, UserOptionsForm, ExperimentForm, ExperimentFlowcellForm

from web.utils import get_run_details, split_flowcell

from django.contrib import messages




def index(request):

    if request.user.is_authenticated:

        return redirect('flowcells')

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
    return render(request, 'web/about_nav.html')


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
    return render(request, 'web/help.html')


@login_required
def flowcells(request):
    flowcells = Flowcell.objects.filter(owner=request.user)
    return render(request, 'web/flowcells.html', context={'flowcells': flowcells})


@login_required
def flowcell_index(request, pk):
    """
    Return the HTML for a flowcells own page after selecting it from the flowcell table.
    :param request: HTTP request object
    :param pk: The primary key of the flowcell
    :type pk: int
    :return: The base HTMl page for looking at a flowcell.
    """
    user = request.user

    flowcell = Flowcell.objects.get(pk=pk)

    if user == flowcell.owner or user.has_perm('view_data', flowcell) or user.has_perm('run_analysis', flowcell):

        return render(request, 'web/flowcell_index.html', context={'flowcell': flowcell})

    else:

        return render(request, 'web/404.html')


@login_required
def flowcell_reads(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    return render(request, 'web/flowcell_reads.html', context={'flowcell': flowcell})


def flowcell_list(request):
    flowcell_list = Flowcell.objects.filter(owner=request.user)
    return render(request, 'web/flowcell_list.html', context={'access_flowcell': flowcell_list})


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

    flowcell_id = int(request.GET.get('flowcell_id', 0))

    draw = int(request.GET.get('draw', 0))

    search_value = request.GET.get('search[value]', '')

    start = int(request.GET.get('start', 0))

    length = int(request.GET.get('length', 10))

    end = start + length
    # Which column is ordering
    order_column = request.GET.get('order[0][column]', '')
    # ascending descending
    order_dir = request.GET.get('order[0][dir]', '')

    run_list = Run.objects.filter(flowcell__id=flowcell_id).filter(flowcell__owner=request.user)

    run_id_list = []

    for run in run_list:

        run_id_list.append(run.id)

    if not search_value == "":

        if search_value[0] == '>':

            reads_temp = FastqRead.objects\
                .filter(run__in=run_id_list)\
                .filter(sequence_length__gt=search_value[1:])
                #.filter(run__flowcell_id=flowcell_id)\
                #.filter(run__flowcell__owner=request.user)\

        elif search_value[0] == '<':

            reads_temp = FastqRead.objects\
                .filter(run__in=run_id_list)\
                .filter(sequence_length__lt=search_value[1:])
                #.filter(run__flowcell_id=flowcell_id) \
                #.filter(run__flowcell__owner=request.user)\

        else:

            reads_temp = FastqRead.objects \
                .filter(run__in=run_id_list)\
                .filter(read_id__contains=search_value)
                #.filter(run__flowcell_id=flowcell_id) \
                #.filter(run__flowcell__owner=request.user) \

    else:

        reads_temp = FastqRead.objects\
            .filter(run__in=run_id_list)
            #.filter(run__flowcell_id=flowcell_id)\
            #.filter(run__flowcell__owner=request.user)\

    if order_column:

        if order_dir == 'desc':

            reads_temp2 = reads_temp.order_by('-{}'.format(query_columns[int(order_column)]))

        else:

            reads_temp2 = reads_temp.order_by('{}'.format(query_columns[int(order_column)]))

    else:
        reads_temp2 = reads_temp

    reads = reads_temp2[start:end]\
        .values('run__runid', 'barcode__name', 'read_id', 'read', 'channel', 'sequence_length')

    records_total = reads_temp.count()

    result = {
        'draw': draw,
        "recordsTotal": records_total,
        "recordsFiltered": records_total,
        "data": list(reads)
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

@login_required
def flowcell_run_stats_download(request, pk):

    flowcell = Flowcell.objects.get(pk=pk)

    minionrunstats_list = MinIONRunStats.objects.filter(run_id__in=flowcell.runs.all())

    runstats_dataframe = pd.DataFrame(list(minionrunstats_list.values()))
    runstats_dataframe = runstats_dataframe.set_index('sample_time')
    runstats_dataframe = runstats_dataframe.drop('created_date', axis=1)
    runstats_dataframe = runstats_dataframe.drop('id', axis=1)
    runstats_dataframe = runstats_dataframe.drop('mean_ratio', axis=1)
    runstats_dataframe = runstats_dataframe.drop('minION_id', axis=1)
    runstats_dataframe = runstats_dataframe.drop('run_id_id', axis=1)
    runstats_dataframe = runstats_dataframe.drop('in_strand', axis=1)
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename={}_live_records.csv'.format(flowcell.name)
    column_order = ['asic_temp',
                    'heat_sink_temp',
                    'voltage_value',
                    'minKNOW_read_count',
                    'event_yield',
                    'above',
                    'adapter',
                    'below',
                    'good_single',
                    'inrange',
                    'multiple',
                    'open_pore',
                    'pore',
                    'zero',
                    'no_pore',
                    'pending_mux_change',
                    'saturated',
                    'strand',
                    'unavailable',
                    'unblocking',
                    'unclassified',
                    'unknown',
                    'minKNOW_histogram_bin_width',
                    'minKNOW_histogram_values']
    runstats_dataframe[column_order].to_csv(response)

    return response


@login_required
def metagenomics_data_download(request, pk):
    """
    Send the Metagenomics data back in CSV format, used by the download button found on the metagenomics tab
    centrifuge/templates/centrifuge/visualisation.html
    :param request:
    :return:
    """

    flowcell = Flowcell.objects.get(pk=pk)
    job_type = JobType.objects.get(name="Metagenomics")
    metagenomics_task = JobMaster.objects.get(flowcell=flowcell, job_type=job_type)
    centrifuge_df = pd.DataFrame(list(CentrifugeOutput.objects.filter(task=metagenomics_task)
                                      .exclude(barcode_name="No").values()))
    centrifuge_df.rename(columns={'classy': 'class'}, inplace=True)
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename={}_centrifuge_output.csv'.format(flowcell.name)
    column_order = ['barcode_name', 'tax_id', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
                    'num_matches', 'proportion_of_classified']
    centrifuge_df[column_order].to_csv(response)
    return response


@login_required
def flowcell_manager(request):
    """
<<<<<<< HEAD
    Returns the flowcell manager HTML
    Parameters
    ----------
    request: django.core.handlers.wsgi.WSGIRequest

    Returns
    -------
    The flowcell manager page HTML template
    """
    # flowcells = Flowcell.objects.filter(owner=request.user)
    return render(request, 'web/flowcell_manager.html', context={'flowcell_manager': flowcell_manager})


def render_messages(request):
    """
    Return the messages html to display messages on the profile tab.
    Parameters
    ----------
    request: django.core.handlers.wsgi.WSGIRequest
        The request object for this AJAX call

    Returns
    -------
    The messages page HTML template
    """
    messages = Message.objects.filter(recipient=request.user).order_by('-created_date')
    flowcells = Message.objects.all().values_list("sender__flowcells__name", flat=True).distinct()


    return render(request, 'web/messages.html', context={"messages": messages, "flowcells": flowcells})


@login_required
def flowcell_manager_runs(request, pk):
    """
    This view shows the a list of runs of a specific flowcell
    in the flowcell maintenance second page
    """

    flowcell = Flowcell.objects.get(pk=pk)
    run_list = get_run_details(pk)
    flowcell_list = Flowcell.objects.filter(owner=request.user).exclude(id=flowcell.id)

    return render(request,
                  'web/flowcell_manager_runs.html',
                  context={'run_list': run_list,
                           'flowcell_list': flowcell_list,
                           'flowcell': flowcell})

@login_required
def flowcell_manager_runs_split(request, pk):
    """
    This view calls the split_flowcell function and
    shows a status page for the user
    """

    flowcell_id = request.POST.get('flowcell_id', None)
    run_id = request.POST.get('run_id', None)
    new_flowcell_name = request.POST.get('new_flowcell_name', None)
    new_or_existing_flowcell = request.POST.get('new_or_existing_flowcell', None)
    existing_flowcell_id = request.POST.get('existing_flowcell_id', None)

    run, from_flowcell, to_flowcell = split_flowcell(new_or_existing_flowcell, flowcell_id, existing_flowcell_id, new_flowcell_name, run_id)

    return render(request,
                  'web/flowcell_manager_runs_split.html',
                  context={'run': run,
                           'from_flowcell': from_flowcell,
                           'to_flowcell': to_flowcell})
