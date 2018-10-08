import json

import numpy as np
from django.db.models import Sum
from django.http import HttpResponse, JsonResponse
from rest_framework.decorators import api_view

from alignment.models import PafRoughCov, PafSummaryCov
# from alignment.models import PafSummaryCov
from alignment.serializers import (PafRoughCovChromSerializer,
                                   )
from alignment.util import calculate_coverage_new
from reference.models import ReferenceLine
from . import util


class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, np.integer):
            return int(obj)

        elif isinstance(obj, np.floating):
            return float(obj)

        else:
            return super(MyEncoder, self).default(obj)


def find_bin(start, size_of_each_bin, value):
    return int((value - start) / size_of_each_bin)


@api_view(['GET'])
def rough_coverage_complete_chromosome_run(request, run_id, barcode_id, read_type_id, chromosome_id):
    return paf_alignment_list(request, run_id, barcode_id, read_type_id, chromosome_id, 0, 0)


@api_view(['GET'])
def paf_alignment_list(request, run_id, barcode_id, read_type_id, chromosome_id, start, end):

    min_extreme = request.GET.get('min', '')
    max_extreme = request.GET.get('max', '')

    if min_extreme != '' and max_extreme != '':

        print('Running query with min and max {} {}'.format(min_extreme, max_extreme))

        min_extreme = int(min_extreme)
        max_extreme = int(max_extreme)

        if min_extreme < 0:
            min_extreme = 0

        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(barcode__id=barcode_id) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .filter(p__gte=min_extreme) \
            .filter(p__lte=max_extreme) \
            .order_by('p')

        current_incdel = get_incdel_at_position(run_id, barcode_id, read_type_id, chromosome_id, min_extreme, False)

    else:

        print('Running query without min and max')

        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(barcode__id=barcode_id) \
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

    for item in queryset:
        positions.append(item.p)
        incdels.append(item.i)

    bin_results = util.calculate_coverage(positions, incdels, current_incdel, reference_size, number_of_bins, bin_width, bin_edges, min_extreme, max_extreme)

    result_list = []

    for key in bin_results.keys():
        result_list.append((min_extreme + (key * bin_width), bin_results[key]))

    return HttpResponse(json.dumps(result_list, cls=NumpyEncoder), content_type="application/json")


@api_view(['GET'])
def rough_coverage_complete_chromosome_flowcell(request, task_id, barcode_name, read_type_id, chromosome_id):

    reference_line = ReferenceLine.objects.get(pk=chromosome_id)

    result_list = calculate_coverage_new(
        request.user,
        task_id,
        barcode_name,
        read_type_id,
        chromosome_id,
        0,
        reference_line.chromosome_length
    )

    # return flowcell_paf_alignment_list(request, task_id, barcode_name, read_type_id, chromosome_id, 0, 0)

# def rough_coverage_complete_chromosome_flowcell(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id):
#     return flowcell_paf_alignment_list(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id, 0, 0)
    return JsonResponse(json.dumps(result_list, cls=NumpyEncoder), safe=False)


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

    bin_results = util.calculate_coverage(positions, incdels, current_incdel, reference_size, number_of_bins, bin_width, bin_edges, min_extreme, max_extreme)

    result_list = []

    for key in bin_results.keys():
        result_list.append((min_extreme + (key * bin_width), bin_results[key]))

    return HttpResponse(json.dumps(result_list, cls=NumpyEncoder), content_type="application/json")


@api_view(['GET'])
def paf_alignment_summary(request, pk):#,bc,ch):
    pass

    """

    :param request:
    :param pk:
    :return:
    """
    # if request.method == 'GET':
    #     queryset = PafSummaryCov.objects \
    #         .filter(run__owner=request.user) \
    #         .filter(run__id=pk)
    #         #.filter(run__id=pk,barcode=bc,chromosome=ch)
    #
    #     serializer = PafSummaryCovSerializer(queryset, many=True, context={'request': request})
    #
    #     return JsonResponse(serializer.data)\


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
def paf_alignment_ch_list(request, pk, ch):

    if request.method == 'GET':

        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)\
            .filter(chromosome__id=ch) \
            .order_by('p')

        serializer = PafRoughCovChromSerializer(queryset, many=True, context={'request': request})

        return JsonResponse(serializer.data)


@api_view(['GET'])
def references_used_by_run(request, run_id):

    references = PafRoughCov.objects.\
        filter(run_id=run_id).\
        values('chromosome__id', 'barcode', 'barcode__name', 'reference', 'reference__reference_name', 'chromosome__line_name', 'read_type', 'read_type__name').distinct()

    result = [{'run_id': run_id, 'chromosome_id': r['chromosome__id'], 'chromosome_name': r['chromosome__line_name'], 'barcode_id': r['barcode'], 'barcode_name': r['barcode__name'], 'reference_id': r['reference'], 'reference_name': r['reference__reference_name'], 'read_type_id': r['read_type'], 'read_type_name': r['read_type__name']} for r in references]
    return HttpResponse(json.dumps(list(result)), content_type="application/json")


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
