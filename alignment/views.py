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


def find_bin(start, size_of_each_bin, value):
    return int((value - start) / size_of_each_bin)


@api_view(['GET'])
def rough_coverage_complete_chromosome(request, run_id, barcode_id, read_type_id, chromosome_id):
    return paf_alignment_list(request, run_id, barcode_id, read_type_id, chromosome_id, 0, 0)


@api_view(['GET'])
def paf_alignment_list(request, run_id, barcode_id, read_type_id, chromosome_id, start, end):

    queryset = PafRoughCov.objects \
        .filter(run__owner=request.user) \
        .filter(barcode__id=barcode_id) \
        .filter(chromosome__id=chromosome_id) \
        .filter(read_type__id=read_type_id) \
        .order_by('p')

    start = int(start)
    end = int(end)

    if start == 0 and end == 0:
        start = queryset.first().p
        end = queryset.last().p

    result_list = []
    current_coverage_sum = 0

    for key, item in enumerate(queryset):

        current_coverage_sum = current_coverage_sum + item.i

        if start < item.p < end:
            result_list.append([item.p, current_coverage_sum])

    return HttpResponse(json.dumps(result_list), content_type="application/json")


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
