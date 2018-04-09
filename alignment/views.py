import json

import numpy as np
import scipy.stats
from django.conf import settings
from django.core.paginator import EmptyPage, PageNotAnInteger, Paginator
from django.db.models import Sum
from django.http import HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import (PafRoughCov, PafSummaryCov,
                              PafSummaryCov_transcriptome)
from alignment.serializers import (PafRoughCovChromSerializer,
                                   PafRoughCovChromSerializerCount,
                                   PafSummaryCovSerializer)
from jobs.models import JobMaster
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


def get_incdel_at_position(run_or_flowcell_id, barcode_or_barcodegroup_id, read_type_id, chromosome_id, position, is_flowcell):

    if is_flowcell:
        queryset = PafRoughCov.objects \
            .filter(flowcell__id=run_or_flowcell_id) \
            .filter(barcodegroup__id=barcode_or_barcodegroup_id) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .filter(p__lte=position) \
            .aggregate(incdel=Sum('i'))

    else:
        queryset = PafRoughCov.objects \
            .filter(run__id=run_id) \
            .filter(barcode__id=barcode_id) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .filter(p__lte=position) \
            .aggregate(incdel=Sum('i'))

    return queryset['incdel']


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
def rough_coverage_complete_chromosome_flowcell(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id):
    return flowcell_paf_alignment_list(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id, 0, 0)


@api_view(['GET'])
def flowcell_paf_alignment_list(request, flowcell_id, barcodegroup_id, read_type_id, chromosome_id, start, end):

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
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafSummaryCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = PafSummaryCovSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)\


@api_view(['GET'])
def flowcellpaf_alignment_summary(request, pk):#,bc,ch):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafSummaryCov.objects \
            .filter(flowcell_id=pk)
            #.filter(run__id=pk,barcode=bc,chromosome=ch)

        serializer = PafSummaryCovSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)\


@api_view(['GET'])
def paf_alignment_transcriptome_summary(request, pk):#,bc,ch):
    """
    CUrrent API behaviour asssumes starting from page1.
    Request for page 0 returns last page.
    This could be a source of confusion and we should resolve.
    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafSummaryCov_transcriptome.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk) \
            .order_by('-read_count')
            #.filter(run__id=pk,barcode=bc,chromosome=ch)


        paginator = Paginator(queryset, settings.PAGINATION_PAGE_SIZE)

        page = request.GET.get('page')

        try:
            result = paginator.page(page)
        except PageNotAnInteger:
            result = paginator.page(1)
        except EmptyPage:
            result = paginator.page(paginator.num_pages)

        serializer = PafSummaryCovSerializer(result, many=True, context={'request': request})

        result2 = {}
        result2['number_pages'] = paginator.num_pages
        result2['data'] = serializer.data
        return HttpResponse(json.dumps(result2), content_type="application/json")

        #serializer = PafSummaryCovSerializer(queryset, many=True, context={'request': request})

        #return Response(serializer.data)


@api_view(['GET'])
def paf_test(request, pk,bc,ch,ty,po,ln):
    """

    :param request:
    :param pk:
    :return:
    """
    if request.method == 'GET':
        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk) \
            .filter(barcode__id=bc) \
            .filter(chromosome__id=ch) \
            .filter(read_type__id=ty) \
            .filter(p__lte=po).values('reference').annotate(sumcount=Sum('i'))
        # .filter(run__id=pk,barcode=bc,chromosome=ch)

        #serializer = PafRoughCovSerializer(queryset, many=True, context={'request': request})

        #queryset = PafRoughCov.objects.filter(run__owner=request.user).filter(run__id=pk).filter(p__lte=1000000).values('reference').annotate(sumcount = Sum('i')).order_by('sumcount')

        serializer1 = PafRoughCovChromSerializerCount(queryset, many=True, context={'request': request})

        queryset2 = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk) \
            .filter(barcode__id=bc) \
            .filter(chromosome__id=ch) \
            .filter(read_type__id=ty) \
            .filter(p__gt=po)\
            .filter(p__lt=(po+ln))\
            .order_by('p')

        serializer2 = PafRoughCovChromSerializer(queryset2, many=True, context={'request': request})

        Serializer_list = [serializer1.data, serializer2.data]

        content = {
            'status': 1,
            'responseCode': status.HTTP_200_OK,
            'data': Serializer_list,

        }

        return Response(content)


@api_view(['GET'])
def paf_alignment_ch_list(request, pk, ch):

    if request.method == 'GET':

        queryset = PafRoughCov.objects \
            .filter(run__owner=request.user) \
            .filter(run__id=pk)\
            .filter(chromosome__id=ch) \
            .order_by('p')

        serializer = PafRoughCovChromSerializer(queryset, many=True, context={'request': request})

        return Response(serializer.data)


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
        values('chromosome__id', 'barcode__barcodegroup', 'barcode__name', 'reference', 'reference__reference_name', 'chromosome__line_name', 'read_type', 'read_type__name').distinct()

    result = [{'flowcell_id': flowcell_id, 'chromosome_id': r['chromosome__id'], 'chromosome_name': r['chromosome__line_name'], 'barcode_id': r['barcode__barcodegroup'], 'barcode_name': r['barcode__name'], 'reference_id': r['reference'], 'reference_name': r['reference__reference_name'], 'read_type_id': r['read_type'], 'read_type_name': r['read_type__name']} for r in references]
    return HttpResponse(json.dumps(list(result)), content_type="application/json")
