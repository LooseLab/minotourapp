import json

import pandas as pd
from celery import task
from django.core.cache import cache

from devices.models import Flowcell
from reads.models import HistogramSummary


@task
def flowcell_histogram_summary():

    flowcells = Flowcell.objects.all()

    for flowcell in flowcells:

        #flowcell = Flowcell.objects.get(pk=pk)

        run_list = flowcell.runs.all()

        qs = HistogramSummary.objects\
           .filter(run__in=run_list)\
           .order_by('read_type', 'bin_index')
           #.filter(run__owner=request.user)\

        if len(qs) > 0:
            df = pd.DataFrame.from_records(qs.values('barcode__name', 'read_type__name', 'is_pass', 'read_count', 'read_length', 'bin_index'))

            gb = df.groupby(['barcode__name', 'read_type__name', 'is_pass', 'bin_index']).agg({'read_count': ['sum'], 'read_length': ['sum']})
            # # ## new_index should get the maximum value in the data to rebuild the histogram dataset
            new_index = list(range(0, gb.index.values.max()[3] + 1))
            ## Here we unstack the multiindex datafram so we can reindex the bin_index.
            ## We do this using the new_index and we fill in the blanks with 0 values.
            ## We then restack the data to rebuild the multindex.
            ## Finally we reorder the indexes to return the order to the original
            gb = gb.unstack(['barcode__name', 'read_type__name', 'is_pass']).reindex(new_index, fill_value=0).stack(['is_pass', 'read_type__name', 'barcode__name']).reorder_levels(
                 ['barcode__name', 'read_type__name', 'is_pass', 'bin_index']).sortlevel(level=0)

            payload = gb.reset_index().apply(lambda row: (row['barcode__name'][0], '{} {} {}'.format(row['barcode__name'][0], row['read_type__name'][0], 'Pass' if row['is_pass'][0] == True else 'Fail'), row['bin_index'][0] * 900 + 900, row['read_count']['sum'], row['read_length']['sum']), axis=1)

            indexes = gb.reset_index().apply(lambda row: '{} {} {}'.format(row['barcode__name'][0], row['read_type__name'][0], 'Pass' if row['is_pass'][0] == True else 'Fail'), axis=1)

            #categories = list(range(0,100000))

            key = "histogram-summary-flowcell-{}".format(flowcell.id)
            cache.set(key, gb.to_json())
            print("Caching {}".format(key))
            # return Response({'data': payload, 'indexes': set(indexes)}) #, 'categories': categories})
            #return Response({'data': gb.to_records()}) #, 'categories': categories})
