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


@api_view(['GET'])
def paf_alignment_list(request, run_id, barcode_id, read_type_id, chromosome_id):

    queryset = PafRoughCov.objects \
        .filter(run__owner=request.user) \
        .filter(barcode__id=barcode_id) \
        .filter(chromosome__id=chromosome_id) \
        .filter(read_type__id=read_type_id) \
        .order_by('p')

    result_list = []
    position_list = []
    coverage_list = []
    current_coverage_sum = 0

    for key, item in enumerate(queryset):

        current_coverage_sum = current_coverage_sum + item.i

        position_list.append(item.p)

        coverage_list.append(current_coverage_sum)

        result_list.append([item.p, current_coverage_sum])

    return HttpResponse(json.dumps(result_list), content_type="application/json")

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

    #jobs = JobMaster.objects.filter(run_id=run_id).filter(job_type__name='Minimap2')

    #references = [job.reference for job in jobs]

    #reference_lines = ReferenceLine.objects.filter(reference__in=references)

    #result = [(reference_line.id, reference_line.reference.reference_name, reference_line.line_name) for reference_line in reference_lines]

    #for job in jobs:
    #    for reference_line in job.reference.referencelines:
    #        print(reference_line.name)
    #        names.append(reference_line.name)

    references = PafRoughCov.objects.\
        filter(run_id=run_id).\
        values('chromosome__id', 'barcode', 'barcode__name', 'reference', 'reference__reference_name', 'chromosome__line_name').distinct()

    print(references[0].keys())

    result = [{'chromosome__id': r['chromosome__id'], 'chromosome__name': r['chromosome__line_name'],'barcode': r['barcode'], 'barcode__name': r['barcode__name'], 'reference': r['reference'], 'reference_name': r['reference__reference_name']} for r in references]
    return HttpResponse(json.dumps(list(result)), content_type="application/json")
