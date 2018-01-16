import json

from django.conf import settings
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db.models import Sum
from django.http import HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from alignment.models import PafRoughCov
from alignment.models import PafSummaryCov, PafSummaryCov_transcriptome
from alignment.serializers import PafRoughCovChromSerializer
from alignment.serializers import PafRoughCovChromSerializerCount
from alignment.serializers import PafSummaryCovSerializer
from reads.models import JobMaster
from reference.models import ReferenceLine

import scipy.stats
import numpy as np


class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, np.integer):
            return int(obj)

        elif isinstance(obj, np.floating):
            return float(obj)

        else:
            return super(MyEncoder, self).default(obj)


def get_incdel_at_position(run_id, barcode_id, read_type_id, chromosome_id, position):

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
def rough_coverage_complete_chromosome(request, run_id, barcode_id, read_type_id, chromosome_id):
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

        current_incdel = get_incdel_at_position(run_id, barcode_id, read_type_id, chromosome_id, min_extreme)

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

    # print(positions)
    # print(incdels)

    # q = np.array([
    #     [3, 10, 15, 23, 70, 99],  # positions
    #     [1, 1, -1, 7, -1, -3]  # incdels
    # ])

    q = np.array([
        positions,
        incdels
    ])

    #
    # this array indicates the bin of each element in the array q[0] (positions)
    # it has the same length of q[0]
    #
    bin_number = np.array([int((x - min_extreme) / bin_width) for x in q[0]])
    # print(bin_number)

    #
    # keep track of current incdel value
    #
    #current_incdel = 0

    bin_results = {}

    #
    # for each bin
    #
    for b in range(number_of_bins):
        positions = q[0][bin_number == b]
        incdels = q[1][bin_number == b]

        begin_bin = min_extreme + (bin_width * b)
        end_bin = min_extreme + (bin_width * (b + 1))

        #
        # contribution is the amount o value that each individual in a bin_number
        # adds to the average and it is the individual's incdel * individual's length
        #
        current_contribution = 0.0


        number_of_individuals = sum(bin_number == b)

        last_individual_index = number_of_individuals - 1

        for i in range(number_of_individuals):

            #
            # calculate the contribution of the last individual of previous bin
            # to the current bin. it happens when a length of an individual crosses
            # the limit of a bin (very very common)
            #
            if i == 0:
                length = positions[0] - begin_bin

                current_contribution = current_contribution + current_incdel * length

                # print(
                # 'bin: {}, i: -, position: {}, incdel: -, current_incdel: {}, length: {}, contribution: {}'.format(b,
                #                                                                                                   begin_bin,
                #                                                                                                   current_incdel,
                #                                                                                                   length,
                #                                                                                                   current_incdel * length))

            #
            # if it is the last individual, length is calculated from its position to
            # the end of the bin
            #
            if i == last_individual_index:
                length = end_bin - positions[i]
            #
            # if it is not the last individual, just calculates its distance until the next
            #
            else:
                length = positions[i + 1] - positions[i]

            #
            # increment incdel
            #
            current_incdel = current_incdel + incdels[i]

            #
            # increment contribution: incdel * length
            #
            current_contribution = current_contribution + current_incdel * length

            # print(
            # 'bin: {}, i: {}, position: {}, incdel: {}, current_incdel: {}, length: {}, contribution: {}'.format(b, i,
            #                                                                                                     positions[
            #                                                                                                         i],
            #                                                                                                     incdels[
            #                                                                                                         i],
            #                                                                                                     current_incdel,
            #                                                                                                     length,
            #                                                                                                     current_incdel * length))

        #
        # if there are individuals, bin average is the sum of the contributions divided by bin width
        #
        if number_of_individuals > 0:
            bin_average = current_contribution / bin_width
            # print('bin average: {}, number of individuals: {}, current contribution: {}'.format(bin_average, number_of_individuals, current_contribution))

        #
        # if bin is empty, then bin average is the value of current incdel (last individual) divided by bin width
        #
        else:
            bin_average = current_incdel
            # print('bin average: {}'.format(bin_average))


        bin_results[b] = bin_average

        # print('bin: {}, average: {:.2f}'.format(b, bin_average))

    # print('bin_width: {}'.format(bin_width))
    # print('bin averages: ')
    # print(bin_results)

    result_list = []

    for key in bin_results.keys():
        result_list.append((min_extreme + (key * bin_width), bin_results[key]))

    return HttpResponse(json.dumps(result_list, cls=NumpyEncoder), content_type="application/json")


@api_view(['GET'])
def paf_alignment_list_old(request, run_id, barcode_id, read_type_id, chromosome_id, start, end):

    NUMBER_OF_BINS = 20

    queryset = PafRoughCov.objects \
        .filter(run__owner=request.user) \
        .filter(barcode__id=barcode_id) \
        .filter(chromosome__id=chromosome_id) \
        .filter(read_type__id=read_type_id) \
        .order_by('p')

    start = int(start)
    end = int(end)

    if start != 0 and end != 0:
        # start = int(start - start * 0.1)
        # end = int(end + end * 0.1)
        number_of_bases_in_a_bin = int((end - start) / NUMBER_OF_BINS)
    else:
        start = queryset.first().p
        end = queryset.last().p
        number_of_bases_in_a_bin = int((end - start) / NUMBER_OF_BINS)

    # print('number of bases in a bin: {} {} {}'.format(number_of_bases_in_a_bin, start, end))

    result_bins = {}

    for i in range(NUMBER_OF_BINS + 1):
        result_bins[i] = {
            'number_of_points': 0,
            'sum_coverage': 0,
            'p': start + (i * number_of_bases_in_a_bin)
        }

    result_list = []
    position_list = []
    coverage_list = []
    current_coverage_sum = 0

    result_per_bin = {}
    result_per_bin2 = {}
    # is_a_new_bin = True
    # current_bin_start_point = None

    for key, item in enumerate(queryset):

        if key == 0:
            # bins = {}
            # current_bin_start_point = item.p
            serie_start_point = item.p

        current_coverage_sum = current_coverage_sum + item.i

        position_list.append(item.p)

        coverage_list.append(current_coverage_sum)

        result_list.append([item.p, current_coverage_sum])

        if item.p > end:
            break # do not process after the end of the screen TODO filter query?

        if item.p >= start:
            bin_index2 = find_bin(start, number_of_bases_in_a_bin, item.p)

            # print('bin index 2: {} {}'.format(bin_index2, current_coverage_sum))

            result_bins[bin_index2]['number_of_points'] =+ 1
            result_bins[bin_index2]['sum_coverage'] =+ current_coverage_sum
            # print(result_bins[bin_index2])

    result_per_bin[queryset.last().p] = current_coverage_sum
    result_per_bin2[queryset.last().p] = current_coverage_sum

    result_list2 = []
    for i in range(NUMBER_OF_BINS):

        p = result_bins[i]['p']

        if result_bins[i]['number_of_points'] >= 1:
            coverage = result_bins[i]['sum_coverage'] / result_bins[i]['number_of_points']

        else:
            coverage = 0

        result_list2.append([p, coverage])

    result = {}
#    result['data_original'] = result_list
    result['data_simplified'] = result_list2


    #return HttpResponse(json.dumps(result_list), content_type="application/json")
    return HttpResponse(json.dumps(result_list2), content_type="application/json")



@api_view(['GET'])
def flowcellpaf_alignment_list(request, run_id, barcodegroup_id, read_type_id, chromosome_id, start, end):

    NUMBER_OF_BINS = 200



    queryset = PafRoughCov.objects \
        .filter(flowcell__owner=request.user) \
        .filter(flowcell_id=run_id)\
        .filter(barcodegroup__id=barcodegroup_id) \
        .filter(chromosome__id=chromosome_id) \
        .filter(read_type__id=read_type_id) \
        .order_by('p')

    start = int(start)
    end = int(end)

    if start != 0 and end != 0:
        # start = int(start - start * 0.1)
        # end = int(end + end * 0.1)
        number_of_bases_in_a_bin = int((end - start) / NUMBER_OF_BINS)
    else:
        start = queryset.first().p
        end = queryset.last().p
        number_of_bases_in_a_bin = int((end - start) / NUMBER_OF_BINS)

    # print('number of bases in a bin: {} {} {}'.format(number_of_bases_in_a_bin, start, end))

    result_bins = {}

    for i in range(NUMBER_OF_BINS + 1):
        result_bins[i] = {
            'number_of_points': 0,
            'sum_coverage': 0,
            'p': start + (i * number_of_bases_in_a_bin)
        }

    result_list = []
    position_list = []
    coverage_list = []
    current_coverage_sum = 0

    result_per_bin = {}
    result_per_bin2 = {}
    # is_a_new_bin = True
    # current_bin_start_point = None

    for key, item in enumerate(queryset):

        if key == 0:
            # bins = {}
            # current_bin_start_point = item.p
            serie_start_point = item.p

        current_coverage_sum = current_coverage_sum + item.i

        position_list.append(item.p)

        coverage_list.append(current_coverage_sum)

        result_list.append([item.p, current_coverage_sum])

        if item.p > end:
            break # do not process after the end of the screen TODO filter query?

        if item.p >= start:
            bin_index2 = find_bin(start, number_of_bases_in_a_bin, item.p)

            # print('bin index 2: {} {}'.format(bin_index2, current_coverage_sum))

            result_bins[bin_index2]['number_of_points'] =+ 1
            result_bins[bin_index2]['sum_coverage'] =+ current_coverage_sum
            # print(result_bins[bin_index2])

    result_per_bin[queryset.last().p] = current_coverage_sum
    result_per_bin2[queryset.last().p] = current_coverage_sum

    result_list2 = []
    for i in range(NUMBER_OF_BINS):

        p = result_bins[i]['p']

        if result_bins[i]['number_of_points'] >= 1:
            coverage = result_bins[i]['sum_coverage'] / result_bins[i]['number_of_points']

        else:
            coverage = 0

        result_list2.append([p, coverage])

    result = {}
#    result['data_original'] = result_list
    result['data_simplified'] = result_list2


    #return HttpResponse(json.dumps(result_list), content_type="application/json")
    return HttpResponse(json.dumps(result), content_type="application/json")



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
