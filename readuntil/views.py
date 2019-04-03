from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from readuntil.serializers import RejectedFastqReadSerializer
import numpy as np
import pickle as picklerick
import pandas as pd


def read_pickle(path_list):
    """
    Read the pickled data
    :param path_list: List of paths to pickle files containing data
    :return:
    """
    listy_list = []
    for path in path_list:
        with open(path, "rb") as fh:
            benefit_dict = picklerick.load(fh)
            # Get the chromosome, which are the keys
            chromosomes = list(benefit_dict.keys())
            # See what chromosomes we have
            print(chromosomes)
            # Loop through them to get the data for each chromosome
            for chromsome in chromosomes:
                a = benefit_dict[chromsome]
                if type(a) is dict:
                    print("dict")
                else:
                    benefit_dict[chromsome] = {"local_benefit": a}

                if listy_list:
                    listy_list[0].update(benefit_dict[chromsome])
                else:
                    listy_list.append(benefit_dict[chromsome])
    return listy_list


def remove_duplicate_sequences(array, master=False):
    """
    Remove stretches of the reference with the same values
    :param array: The array with a value for each reference location
    :param master: Whether this data is for the master chart or not
    :return: a tuple of a list with an x,y coordinate for plotting, the max value on the x axis,
    and the max value on the y axis
    """
    # max on the y value
    ymax = array.max()

    df = pd.DataFrame(array)

    if pd.api.types.is_bool_dtype(df[0]):
        df[0] = df[0].astype(int)
    # max on the x axis, i.e the final base of the reference
    xmax = df.index.values.max()
    # Off set the sequence by one
    df[1] = df[0].shift()
    # Set the new extra value to 0 instead of NaN
    df[1].iloc[0] = 0
    # Remove any values that are the same
    df = df[df[0] != df[1]]
    # Make a list of tuple with the x axis (The index value) and the y axis value (The array value)
    x_y_values = list(zip(df.index.values, df[0]))

    if len(x_y_values) > 50000 and master:
        print("Stormy storm storm storm")
        x_y_values = sampler(x_y_values)

    return x_y_values, xmax, ymax


def get_y_axis_max(array):
    """
    Get the maximum value in the array
    :param array: The array of data
    :type array: np.nparray
    :return:
    """
    return array.max()


def sampler(array):
    """
    If the array is too large > 50000 points, take a sample using steps
    :param array: Numpy array
    :return: a smaller sampled numpy array
    """
    step_size = int(len(array) / 50000)
    sampled_array = array[0::step_size]
    return sampled_array


@api_view(['POST'])
def rejectedfastqread_list(request):

    serializer = RejectedFastqReadSerializer(
        data=request.data,
        context={'request': request}
    )

    if serializer.is_valid():

        serializer.save(owner=request.user)
        return Response(serializer.data, status=status.HTTP_201_CREATED)

    return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


@api_view(["GET"])
def get_benefit_data(request):
    """
    Return data for mock scatter plots
    :param request:
    :return:
    """
    file_list = ["/home/rory/Downloads/rollingdict.p", "/home/rory/Downloads/referencedict.p",
                 "/home/rory/Downloads/benefitdict.p"]

    # file_list = ["/home/rory/Downloads/referencedict.p", "/home/rory/Downloads/benefitdict.p"]

    # file_list = ["/home/rory/Downloads/benefitdict.p"]

    data_array = read_pickle(file_list)

    coverage = data_array[0]["coverage"]

    x_y_cov, xmax_cov, ymax_cov = remove_duplicate_sequences(coverage, True)

    y_max_match = get_y_axis_max(data_array[0]["match"])

    y_max_mismatch = get_y_axis_max(data_array[0]["mismatch"])

    y_max_localben = get_y_axis_max(data_array[0]["local_benefit"])

    y_max_forward_roll = get_y_axis_max(data_array[0]["Forward"])

    y_max_reverse_roll = get_y_axis_max(data_array[0]["Reverse"])

    y_max_fwd_mask = get_y_axis_max(data_array[0]["ForMask"])

    y_max_rev_mask = get_y_axis_max(data_array[0]["RevMask"])

    data_dict = {"coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
                 "match": {"ymax": y_max_match},
                 "mismatch": {"ymax": y_max_mismatch},
                 "localBen": {"ymax": y_max_localben},
                 "forwardRoll": {"ymax": y_max_forward_roll},
                 "reverseRoll": {"ymax": y_max_reverse_roll},
                 "fwdMask": {"ymax": y_max_fwd_mask},
                 "revMask": {"ymax": y_max_rev_mask}}

    return Response(data_dict)


@api_view(["GET"])
def get_benefit_data_complete(request):
    """
    The complete benefit data for the selected regions of the master chart
    :param request: Contains the min and max x axis values in the body
    :return: a dictionary of the data from the pickle for all the detail charts
    """
    mini = int(request.GET.get("min"))
    maxi = int(request.GET.get("max"))

    file_list = ["/home/rory/Downloads/rollingdict.p", "/home/rory/Downloads/referencedict.p",
                 "/home/rory/Downloads/benefitdict.p"]

    # file_list = ["/home/rory/Downloads/referencedict.p", "/home/rory/Downloads/benefitdict.p"]

    # file_list = ["/home/rory/Downloads/benefitdict.p"]

    data_array = read_pickle(file_list)

    coverage = data_array[0]["coverage"][mini:maxi]

    matches = data_array[0]["match"][mini:maxi]

    mismatches = data_array[0]["mismatch"][mini:maxi]

    local_benefit = data_array[0]["local_benefit"][mini:maxi]
    # Forward strand rolling benefit
    forward_roll = data_array[0]["Forward"][mini:maxi]
    # Reverse strand rolling benefit
    reverse_roll = data_array[0]["Reverse"][mini:maxi]

    fwd_mask = data_array[0]["ForMask"][mini:maxi]

    rev_mask = data_array[0]["RevMask"][mini:maxi]

    x_y_cov, xmax_cov, ymax_cov = remove_duplicate_sequences(coverage)

    x_y_match, xmax_match, ymax_match = remove_duplicate_sequences(matches)

    x_y_mismatch, xmax_mismatch, ymax_mismatch = remove_duplicate_sequences(mismatches)

    x_y_local_ben, xmax_local_ben, ymax_local_ben = remove_duplicate_sequences(local_benefit)

    x_y_fwd_rolling_ben, xmax_fwd_rolling_ben, ymax_fwd_rolling_ben = remove_duplicate_sequences(forward_roll)

    x_y_rev_rolling_ben, xmax_rev_rolling_ben, ymax_rev_rolling_ben = remove_duplicate_sequences(reverse_roll)

    x_y_fwd_mask, xmax_fwd_mask, ymax_fwd_mask = remove_duplicate_sequences(fwd_mask)

    x_y_rev_mask, xmax_rev_mask, ymax_rev_mask = remove_duplicate_sequences(rev_mask)

    data_dict = {"coverage": {"xmax": xmax_cov, "ymax": ymax_cov, "data": x_y_cov},
                 "match": {"xmax": xmax_match, "ymax": ymax_match, "data": x_y_match},
                 "mismatch": {"xmax": xmax_mismatch, "ymax": ymax_mismatch, "data": x_y_mismatch},
                 "localBenefit": {"xmax": xmax_local_ben, "ymax": ymax_local_ben, "data": x_y_local_ben},
                 "rollingBenefitFwd": {"xmax": xmax_fwd_rolling_ben, "ymax": ymax_fwd_rolling_ben,
                                       "data": x_y_fwd_rolling_ben},
                 "rollingBenefitRev": {"xmax": xmax_rev_rolling_ben, "ymax": ymax_fwd_rolling_ben,
                                       "data": x_y_rev_rolling_ben},
                 "forwardMask": {"xmax": xmax_fwd_mask, "ymax": ymax_fwd_mask,
                                 "data": x_y_fwd_mask},
                 "reverseMask": {"xmax": xmax_rev_mask, "ymax": ymax_rev_mask,
                                 "data": x_y_rev_mask},
                 }

    return Response(data_dict)


