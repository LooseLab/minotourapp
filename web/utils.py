import collections
import datetime
import math
import os
import sys
from functools import lru_cache

import numpy as np

from minotourapp.utils import get_env_variable
from reads.models import Flowcell, MinionRunInfo


def file_folder_exists(filepath):
    """
    Check if a provided file or folder exists.

    Parameters
    ----------
    filepath: str
        Absolute path to file or folder to be tested.

    Returns
    -------
    bool

    """
    return os.path.exists(filepath)


def is_executable(filepath):
    """
    Check whether the provided file path is executable.
    Parameters
    ----------
    filepath: str

    Returns
    -------
    bool

    """
    return os.access(filepath, os.X_OK)


def validate_envs():
    """
    Validate environment variables for obvious errors.
    Returns
    -------
    None
    TODO finish this
    """
    # Artic checks
    if " " in get_env_variable("MT_ARTIC_RESULTS_DIR"):
        raise ValueError("Please remove whitespace from the Artic results directory variable name.")
    if not file_folder_exists(get_env_variable("MT_ARTIC_RESULTS_DIR")):
        raise ValueError("The Artic results directory specified in the environment variables file does not exist.")

    # Minimap2 checks
    if not is_executable(get_env_variable("MT_MINIMAP2")):
        raise ValueError("The Minimap2 executable is either not there, or not executable. Please double check.")


def update_last_activity_time(flowcell):
    """
    Update the last activity date on the flowcell
    :param flowcell: The flowcell to be updated
    :return:
    """
    flowcell.last_activity_date = datetime.datetime.now(datetime.timezone.utc)
    flowcell.save()


def _rstrip_tuple(s, chars):
    """
    Return a tuple of the initial string and the chars removed from the right
    """
    h = s.rstrip(chars)
    return h, s[len(h) :]


def int_string(i):
    """
    Converts integers from string, but leaves strings alone
    https://www.youtube.com/watch?v=D7ImcrILvEo

    Parameters
    ----------
    i : str
        A string that maybe an integer

    Returns
    -------
    int or str
    """
    try:
        return int(i)
    except ValueError:
        return i


def md_generator(s):
    """Given an MD string yields elements in reverse order

    Parameters
    ----------
    s : str
        A "MD:Z" string from a mapping

    Returns
    -------
    Generator that partions the MD tag in reverse order
    """
    while s:
        s, e = _rstrip_tuple(s, "0123456789")
        if e == "":
            s, e = _rstrip_tuple(s, "ATCGN^")
        yield int_string(e)


def parse_MD(s, match=1, mismatch=0):
    """
    Parse an MD string into a list of integers
    Parameters
    ----------
    s : str
        A "MD:Z" string from a mapping

    Returns
    -------
    list
        Returns an MD string as integers, with 1 for a match and 0 for a mismatch
    """
    x = []
    for i in md_generator(s):
        if isinstance(i, int):
            x.extend([match] * i)
        # elif i[0] == "^":
        #    x.extend([mismatch] * (len(i) - 1))
        else:
            x.extend([mismatch] * len(i.replace("^", "")))
    return x[::-1]


def parse_mdpaf_alex(pafline):
    """
    This function takes a line from a paf file created by minimap2 and returns a number of outputs.
    Parameters
    ----------
    pafline: str
        A line from a paf file.

    Returns
    -------
    mismatcharray:
        A numpy array of the same length as the alignment with a 1 for a mismatch and a 0 for a match to the reference.
    matArr:
        A numpy array inverse to the mutArr - so 1 for a match and 0 for a mismatch
    mapstart:
        The starting coordinate for a match to the reference (note: 0 based)
    maporientation:
        Either + or - for strand.
    reference:
        The string identifying the reference ID.
    referencelength:
        The length of the reference sequence
    readlength:
        The length of the read - note this may be longer than the alignement length.
    """

    bits = pafline.split()
    reference = bits[5]
    referencelength = bits[6]
    readlength = bits[1]
    sub = "MD:Z:"
    mapstart = int(bits[7])
    mapend = int(bits[8])
    maporientation = bits[4]
    mdstr = ""

    for s in filter(lambda x: sub in x, bits):
        mdstr = s.split(":")[-1]

    matcharray = np.fromiter(parse_MD(mdstr), np.int)
    mismatcharray = np.ones(matcharray.shape) - matcharray

    return (
        mismatcharray,
        matcharray,
        mapstart,
        mapend,
        maporientation,
        reference,
        referencelength,
        readlength,
    )


def add_chromosome(
    referencedict,
    rollingdict,
    benefitdict,
    reference,
    referencelength,
    mean,
    error,
    prior_diff,
):

    reference_length = int(referencelength)

    referencedict[reference] = dict()
    referencedict[reference]["length"] = reference_length
    referencedict[reference]["mismatch"] = np.zeros(reference_length)
    referencedict[reference]["match"] = np.zeros(reference_length)
    referencedict[reference]["coverage"] = np.zeros(reference_length)

    benefitdict[reference] = getbenefit(
        referencedict[reference]["match"],
        referencedict[reference]["mismatch"],
        error=error,
        prior_diff=prior_diff,
    )

    a = benefitdict[reference]

    rollingdict[reference] = dict()
    rollingdict[reference]["Forward"] = rolling_sum([a], n=int(np.floor(mean)))[
        0
    ]  # what does it do?
    rollingdict[reference]["Reverse"] = rolling_sum([a[::-1]], n=int(np.floor(mean)))[
        0
    ][::-1]
    rollingdict[reference]["ForMean"] = rollingdict[reference]["Forward"].mean()
    rollingdict[reference]["RevMean"] = rollingdict[reference]["Reverse"].mean()
    rollingdict[reference]["ForMask"] = check_mask(rollingdict[reference]["Forward"], 0)
    rollingdict[reference]["RevMask"] = check_mask(rollingdict[reference]["Reverse"], 0)

    return referencedict, rollingdict


### Helper functions to calculate expected benefit, rolling sum and the check mask.


def getbenefit(matchholdarray, mismatchholdarray, error=0.10, prior_diff=0.0001):
    """
    This function takes in two arrays - a count of matches and mismatches at each position and calculates the benefit at each position.

    Parameters:
    -----------
    matchholdarray:
        An array with counts of matches at every position.
    mismatchholdarray:
        An array with counts of mismatches at every position.

    Returns:
    benefitarray:
        An array with the benefit for each position.
    """
    benefit = list()
    for i in range(len(matchholdarray)):
        ### This block of code will artificially inflate the benefit of sequencing a region with zero cove
        # if matchholdarray[i]==0 and mismatchholdarray[i]==0:
        #    benefit.append(0.45)
        #    continue
        benefit.append(
            match_prob(
                matchholdarray[i],
                mismatchholdarray[i],
                error=error,
                prior_diff=prior_diff,
            )
        )
        if np.isnan(benefit[-1]):
            print(
                "We've got a nan {} {}".format(matchholdarray[i], mismatchholdarray[i])
            )
    benefitarray = np.array(benefit)
    return benefitarray


def rolling_sum(a, n=100, pad=True):
    """
    This function returns a rolling sum across an array with a window size of n.
    If pad is true, the input array is padded with zeros at the end to enable calculating the benefit at the -1 position of the array.

    Parameters:
    -----------
    a
        The array over which the rolling sum will be calculated.
    n
        The window size to calculate the rolling sum over. Assumed to be the mean read length.
    pad
        Default to True. Enables calculating the rolling sum to the end of the array.

    Returns:
    --------
        An array of the rolling sum of the input.
    """

    if pad:
        # print ("Padding with n={}".format(n))
        a[0] = np.concatenate((a[0], np.zeros(n - 1)))
    ret = np.cumsum(a, axis=1, dtype=float)
    ret[:, n:] = ret[:, n:] - ret[:, :-n]
    return ret[:, n - 1 :]


def check_mask(benefitarray, mean):
    """
    Function to generate a mask of 1 and 0 for sequencing or rejecting a read.
    At present works by asking if the value of the array is greater than the mean.
    This is problematic as it ends up rejecting a lot of stuff that we want to sequence when coverage is low.
    So - perhaps only update those positions for which there is coverage?

    :param benefitarray: Array of values
    :param mean: mean value of benefit array
    :return: Array of 1 and 0 where 1 will mean sequence and 0 will mean reject - probably represent as bools?
    """
    ### Note - we might want change this to median, not mean.
    return benefitarray >= mean


## This block of functions is used to calculate the expected benefit at a given position given matches, mismatches, prior_diff and error.
## The use of the lru_cache function effectively allows us to cache the calculations so given the same input the function looks up the answer rather than recalculating it.
## I believe there are errors in the code below!


@lru_cache(maxsize=None)
def match_prob(match, mismatch, prior_diff=0.0001, error=0.01):
    non_error = 1 - error
    match_prob = (
        return_posterior_ref(match, mismatch, prior_diff=prior_diff, error=error)
        * non_error
    ) + (
        return_posterior_alt(match, mismatch, prior_diff=prior_diff, error=error)
        * error
    )

    mismatch_prob = (
        return_posterior_alt(match, mismatch, prior_diff=prior_diff, error=error)
        * non_error
    ) + (
        return_posterior_ref(match, mismatch, prior_diff=prior_diff, error=error)
        * error
    )

    matchp1 = match + 1

    mismatchp1 = mismatch + 1

    before = (
        return_posterior_ref(match, mismatch, prior_diff=prior_diff, error=error),
        return_posterior_alt(match, mismatch, prior_diff=prior_diff, error=error),
    )

    after_match = (
        return_posterior_ref(matchp1, mismatch, prior_diff=prior_diff, error=error),
        return_posterior_alt(matchp1, mismatch, prior_diff=prior_diff, error=error),
    )

    after_mismatch = (
        return_posterior_ref(match, mismatchp1, prior_diff=prior_diff, error=error),
        return_posterior_alt(match, mismatchp1, prior_diff=prior_diff, error=error),
    )

    weight = match_prob * kldistance(after_match, before)
    weight += mismatch_prob * kldistance(after_mismatch, before)

    return round(weight, 10)


@lru_cache(maxsize=None)
def return_posterior_ref(match, mismatch, prior_diff=0.001, error=0.01):
    non_error = 1 - error
    base_is_ref = (non_error ** match) * (error ** mismatch) * (1 - prior_diff)
    base_is_alt = (error ** match) * (non_error ** mismatch) * (prior_diff)

    # To prevent division by zero at high coverage
    if (base_is_ref + base_is_alt) == 0:
        ref_alt_prob = sys.float_info.min
    else:
        ref_alt_prob = base_is_ref + base_is_alt

    return base_is_ref / ref_alt_prob


@lru_cache(maxsize=None)
def return_posterior_alt(match, mismatch, prior_diff=0.001, error=0.01):
    non_error = 1 - error
    base_is_ref = (non_error ** match) * (error ** mismatch) * (1 - prior_diff)
    base_is_alt = (error ** match) * (non_error ** mismatch) * (prior_diff)

    # To prevent division by zero at high coverage
    if (base_is_ref + base_is_alt) == 0:
        ref_alt_prob = sys.float_info.min
    else:
        ref_alt_prob = base_is_ref + base_is_alt

    return base_is_alt / ref_alt_prob


@lru_cache(maxsize=None)
def kldistance(after, before):
    #
    # P || Q (Q to P ) is Sum P(i) (log(Pi) - log(Qi))
    #
    count = 0
    for i in range(len(after)):
        # To prevent division by zero at high coverage
        if after[i] == 0:
            aft_i = sys.float_info.min
        else:
            aft_i = after[i]
        if before[i] == 0:
            bef_i = sys.float_info.min
        else:
            bef_i = before[i]
        count = count + (aft_i * (math.log(aft_i) - math.log(bef_i)))

    return count


def RepresentsInt(i):
    """
    Helper function that returns True if `i` is an int.
    """
    try:
        int(i)
        return True
    except ValueError:
        return False


def parse_md_cg_pafline(line):
    """
    This function splits the paf line by \t and returns a named tuple
    """

    PafLine = collections.namedtuple("PafLine", "read_id chromosome")

    columns = line.split("\t")

    read_id = columns[0]

    ref = columns[5]

    chromosome_temp = ref.split("|")
    chromosome = chromosome_temp[1]

    paffile = PafLine(read_id=read_id, chromosome=chromosome)

    return paffile


def get_run_details(run_id):
    # TODO i hate this code, let's rewrite this code and comment it
    flowcell = Flowcell.objects.get(pk=run_id)

    result = {}

    for run in flowcell.runs.all():

        element = {
            "id": None,
            "runid": None,
            "run_start_time": None,
            "first_read": None,
            "last_read": None,
            "minknow_computer_name": None,
            "minion_id": None,
            "asic_id": None,
            "sequencing_kit": None,
            "purpose": None,
            "minknow_version": None,
            "flowcell_type": None,
            "flowcell_id": None,
            "sample_name": None,
            "experiment_name": None,
            "read_count": None,
            "total_read_length": None,
            "max_read_length": None,
            "min_read_length": None,
            "avg_read_length": None,
            "first_read_start_time": None,
            "last_read_start_time": None,
        }

        minion_run_status_list = MinionRunInfo.objects.filter(run_id=run).order_by(
            "minKNOW_start_time"
        )

        if len(minion_run_status_list) > 0:

            minion_run_status = minion_run_status_list[0]

            element["id"] = minion_run_status.run.id
            element["runid"] = minion_run_status.run.runid
            element["run_start_time"] = minion_run_status.minKNOW_start_time
            element["minknow_computer_name"] = minion_run_status.minKNOW_computer
            element["minion_id"] = minion_run_status.minion.minION_name
            element["asic_id"] = minion_run_status.minKNOW_asic_id
            element["sequencing_kit"] = minion_run_status.sequencing_kit
            element["purpose"] = minion_run_status.minKNOW_exp_script_purpose
            element["minknow_version"] = minion_run_status.minKNOW_version
            element["flowcell_type"] = minion_run_status.flowcell_type
            element["flowcell_id"] = minion_run_status.minKNOW_flow_cell_id
            element["sample_name"] = minion_run_status.minKNOW_sample_name
            element["experiment_name"] = minion_run_status.experiment_id

            result[element["runid"]] = element

    result_basecalled_summary = []

    for run in flowcell.runs.all():

        if hasattr(run, "summary"):

            runid = run.summary.run.runid

            if runid in result.keys():

                result[runid]["read_count"] = run.summary.read_count
                result[runid]["total_read_length"] = run.summary.total_read_length
                result[runid]["max_read_length"] = run.summary.max_read_length
                result[runid]["min_read_length"] = run.summary.min_read_length
                result[runid]["avg_read_length"] = run.summary.avg_read_length
                result[runid][
                    "first_read_start_time"
                ] = run.summary.first_read_start_time
                result[runid]["last_read_start_time"] = run.summary.last_read_start_time

            else:

                element = {
                    "runid": None,
                    "run_start_time": None,
                    "first_read": None,
                    "last_read": None,
                    "minknow_computer_name": None,
                    "minion_id": None,
                    "asic_id": None,
                    "sequencing_kit": None,
                    "purpose": None,
                    "minknow_version": None,
                    "flowcell_type": None,
                    "flowcell_id": None,
                    "sample_name": None,
                    "experiment_name": None,
                    "read_count": None,
                    "total_read_length": None,
                    "max_read_length": None,
                    "min_read_length": None,
                    "avg_read_length": None,
                    "first_read_start_time": None,
                    "last_read_start_time": None,
                }

                element["id"] = run.id
                element["runid"] = run.summary.run.runid
                element["read_count"] = run.summary.read_count
                element["total_read_length"] = run.summary.total_read_length
                element["max_read_length"] = run.summary.max_read_length
                element["min_read_length"] = run.summary.min_read_length
                element["avg_read_length"] = run.summary.avg_read_length
                element["first_read_start_time"] = run.summary.first_read_start_time
                element["last_read_start_time"] = run.summary.last_read_start_time

                result[element["runid"]] = element

    return result.values()

