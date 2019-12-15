import numpy as np
from django.db.models import Sum

from alignment.models import PafRoughCov


def calculate_coverage(positions, incdels, current_incdel, reference_size, number_of_bins, bin_width, bin_edges, min_extreme, max_extreme):

    #print(positions)
    #print(incdels)

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

                # print('>>> current_contribution: {}'.format(current_contribution))
                # print('>>> current_incdel: {}'.format(current_incdel))
                # print('>>> length: {}'.format(length))

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
    return bin_results


def get_incdel_at_position(task_id, barcode_name, read_type_id, chromosome_id, position):

    """
    print('>>> get_incdel_at_position')
    print('>>> task_id: {}'.format(task_id))
    print('>>> barcode_name: {}'.format(barcode_name))
    print('>>> read_type_id: {}'.format(read_type_id))
    print('>>> chromosome_id: {}'.format(chromosome_id))
    print('>>> position: {}'.format(position))
    """

    queryset = PafRoughCov.objects \
        .filter(job_master__id=task_id) \
        .filter(barcode_name=barcode_name) \
        .filter(chromosome__id=chromosome_id) \
        .filter(read_type__id=read_type_id) \
        .filter(p__lte=position) \
        .aggregate(incdel=Sum('i'))

    if queryset['incdel']:

        return queryset['incdel']

    else:

        return 0

    return queryset['incdel']


def calculate_coverage_new(user, task_id, barcode_name, read_type_id, chromosome_id, min_extreme, max_extreme):

    if min_extreme != '' and max_extreme != '':

        print('Running query with min and max {} {}'.format(min_extreme, max_extreme))

        min_extreme = int(min_extreme)
        max_extreme = int(max_extreme)

        if min_extreme < 0:
            min_extreme = 0

        queryset = PafRoughCov.objects \
            .filter(flowcell__owner=user) \
            .filter(barcode_name=barcode_name) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .filter(p__gte=min_extreme) \
            .filter(p__lte=max_extreme) \
            .order_by('p') \

        return queryset

        current_incdel = get_incdel_at_position(task_id, barcode_name, read_type_id, chromosome_id, min_extreme)
        # current_incdel = get_incdel_at_position(flowcell_id, barcodegroup_id, read_type_id, chromosome_id, min_extreme, True)

    else:

        print('Running query without min and max')

        queryset = PafRoughCov.objects \
            .filter(flowcell__owner=user) \
            .filter(job_master=task_id) \
            .filter(barcode_name=barcode_name) \
            .filter(chromosome__id=chromosome_id) \
            .filter(read_type__id=read_type_id) \
            .order_by('p') \

        return queryset

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

    #print("queryset length: {}".format(len(queryset)))
    for item in queryset:
        positions.append(item.p)
        incdels.append(item.i)

    bin_results = calculate_coverage(
        positions,
        incdels,
        current_incdel,
        reference_size,
        number_of_bins,
        bin_width,
        bin_edges,
        min_extreme,
        max_extreme
    )

    result_list = []

    for key in bin_results.keys():
        result_list.append((min_extreme + (key * bin_width), bin_results[key]))
