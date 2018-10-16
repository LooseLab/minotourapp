import json

import numpy as np
from django.db.models import Sum
from django.http import HttpResponse, JsonResponse
from rest_framework.decorators import api_view

from alignment.models import PafRoughCov, PafSummaryCov
from alignment.utils import calculate_coverage_new
from . import utils


class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, np.integer):
            return int(obj)

        elif isinstance(obj, np.floating):
            return float(obj)

        else:
            return super(MyEncoder, self).default(obj)


@api_view(['GET'])
def rough_coverage_complete_chromosome_flowcell(request, task_id, barcode_name, read_type_id, chromosome_id):

    queryset = PafRoughCov.objects \
        .filter(job_master__id=task_id) \
        .filter(flowcell__owner=request.user) \
        .filter(barcode_name=barcode_name) \
        .filter(chromosome__id=chromosome_id) \
        .filter(read_type__id=read_type_id) \
        .values('p') \
        .annotate(Sum('i')) \
        .order_by('p')

    result = []

    for record in queryset:

        result.append([record['p'], record['i__sum']])

    return JsonResponse(result, safe=False)


@api_view(['GET'])
def flowcell_paf_alignment_list(request, task_id, barcode_name, read_type_id, chromosome_id, start, end):
    # def flowcell_paf_alignment_list(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id, start, end):

    min_extreme = request.GET.get('min', '')
    max_extreme = request.GET.get('max', '')


    result_list = calculate_coverage_new(task_id, barcode_name, read_type_id, chromosome_id, start, end)

    return JsonResponse(json.dumps(result_list, cls=NumpyEncoder), safe=False)


@api_view(['GET'])
def flowcell_paf_alignment_list_original(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id, start, end):

    min_extreme = request.GET.get('min', '')
    max_extreme = request.GET.get('max', '')

    if min_extreme != '' and max_extreme != '':

        print('Running query with min and max {} {}'.format(min_extreme, max_extreme))

        min_extreme = int(min_extreme)
        max_extreme = int(max_extreme)

        if min_extreme < 0:
            min_extreme = 0

        queryset = PafRoughCov.objects \
            .filter(flowcell__owner=request.user) \
            .filter(barcodegroup__id=barcodegroup_id) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .filter(p__gte=min_extreme) \
            .filter(p__lte=max_extreme) \
            .order_by('p')

        current_incdel = get_incdel_at_position(flowcell_id, barcodegroup_id, read_type_id, chromosome_id, min_extreme, True)

    else:

        print('Running query without min and max')

        queryset = PafRoughCov.objects \
            .filter(flowcell__owner=request.user) \
            .filter(barcodegroup__id=barcodegroup_id) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .order_by('p')

        min_extreme = 0
        max_extreme = queryset[0].chromosome.chromosome_length

        current_incdel = 0

    reference_size = max_extreme - min_extreme

    number_of_bins = 200

    bin_width = reference_size / number_of_bins

    #
    # the size of bin_edges is the number of bins + 1
    #
    bin_edges = np.array([min_extreme + bin_width * i for i in range(number_of_bins + 1)])

    positions = []
    incdels = []

    print("queryset length: {}".format(len(queryset)))
    for item in queryset:
        positions.append(item.p)
        incdels.append(item.i)

    bin_results = utils.calculate_coverage(positions, incdels, current_incdel, reference_size, number_of_bins, bin_width, bin_edges, min_extreme, max_extreme)

    result_list = []

    for key in bin_results.keys():
        result_list.append((min_extreme + (key * bin_width), bin_results[key]))

    return HttpResponse(json.dumps(result_list, cls=NumpyEncoder), content_type="application/json")


@api_view(['GET'])
def flowcell_paf_summary_cov(request, pk):

    queryset = PafSummaryCov.objects.filter(job_master__flowcell_id=pk)

    response = []

    for record in queryset:

        paf_summary_cov = {

            'id': record.job_master.id,
            'barcode_name': record.barcode_name,
            'read_type_name': 'Template',
            'reference_line_name': record.reference_line_name,
            'read_count': record.read_count,
            'total_length': record.total_length,
            'chrom_cover': record.total_length/record.reference_line_length,
            'avg_read_len': record.total_length / record.read_count
        }

        response.append(paf_summary_cov)

    return JsonResponse({'data': response})


@api_view(['GET'])
def flowcellreferences_used_by_run(request, flowcell_id):

    references = PafRoughCov.objects.\
        filter(flowcell_id=flowcell_id).\
        values('job_master__id', 'chromosome__id', 'barcode_name', 'reference', 'reference__name', 'chromosome__line_name', 'read_type', 'read_type__name').distinct()

    result = [
        {
            'task_id': r['job_master__id'],
            'flowcell_id': flowcell_id,
            'chromosome_id': r['chromosome__id'],
            'chromosome_name': r['chromosome__line_name'],
            'barcode_name': r['barcode_name'],
            'reference_id': r['reference'],
            'reference_name': r['reference__name'],
            'read_type_id': r['read_type'],
            'read_type_name': r['read_type__name']
        } for r in references
    ]

    return HttpResponse(json.dumps(list(result)), content_type="application/json")
