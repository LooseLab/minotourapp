import numpy as np


def calculate_coverage(positions, incdels, current_incdel, reference_size, number_of_bins, bin_width, bin_edges, min_extreme, max_extreme):
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
    return bin_results
