import re,sys,os

import time

import numpy as np
import collections
from functools import lru_cache
import math


def _rstrip_tuple(s, chars):
    """
    Return a tuple of the initial string and the chars removed from the right
    """
    h = s.rstrip(chars)
    return h, s[len(h):]


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
        mdstr = s.split(':')[-1]

    matcharray = np.fromiter(parse_MD(mdstr), np.int)
    mismatcharray = np.ones(matcharray.shape) - matcharray

    return mismatcharray, matcharray, mapstart, mapend, maporientation, reference, referencelength, readlength


def add_chromosome(referencedict,rollingdict,benefitdict,reference, referencelength,mean):

    error = 0.10
    prior_diff = 0.0001

    referencedict[reference] = dict()
    referencedict[reference]['length'] = int(referencelength)
    referencedict[reference]["mismatch"] = np.zeros(int(referencelength))
    referencedict[reference]["match"] = np.zeros(int(referencelength))
    referencedict[reference]["coverage"] = np.zeros(int(referencelength))
    benefitdict[reference] = getbenefit(referencedict[reference]["match"],
                                             referencedict[reference]["mismatch"]
                                       ,error=error,prior_diff=prior_diff)
    rollingdict[reference] = dict()
    A = benefitdict[reference]
    rollingdict[reference]["Forward"] = rolling_sum([A], n=int(np.floor(mean)))[0]
    rollingdict[reference]["Reverse"] = rolling_sum([A[::-1]], n=int(np.floor(mean)))[0][::-1]
    rollingdict[reference]["ForMean"] = rollingdict[reference]["Forward"].mean()
    rollingdict[reference]["RevMean"] = rollingdict[reference]["Reverse"].mean()
    rollingdict[reference]["ForMask"] = check_mask(rollingdict[reference]["Forward"],0)
    rollingdict[reference]["RevMask"] = check_mask(rollingdict[reference]["Reverse"],0)
    return referencedict,rollingdict


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
        benefit.append(match_prob(matchholdarray[i], mismatchholdarray[i], error=error, prior_diff=prior_diff))
        if np.isnan(benefit[-1]):
            print ("We've got a nan {} {}".format(matchholdarray[i], mismatchholdarray[i]))
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
    return ret[:, n - 1:]


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
    non_error = (1 - error)
    match_prob = (
                         return_posterior_ref(match, mismatch, prior_diff=prior_diff, error=error) * non_error) + (
                         return_posterior_alt(match, mismatch, prior_diff=prior_diff, error=error) * error)

    mismatch_prob = (
                            return_posterior_alt(match, mismatch, prior_diff=prior_diff, error=error) * non_error) + (
                            return_posterior_ref(match, mismatch, prior_diff=prior_diff, error=error) * error)

    matchp1 = match + 1

    mismatchp1 = mismatch + 1

    before = (return_posterior_ref(match, mismatch, prior_diff=prior_diff, error=error),
              return_posterior_alt(match, mismatch, prior_diff=prior_diff, error=error))

    after_match = (return_posterior_ref(matchp1, mismatch, prior_diff=prior_diff, error=error),
                   return_posterior_alt(matchp1, mismatch, prior_diff=prior_diff, error=error))

    after_mismatch = (return_posterior_ref(match, mismatchp1, prior_diff=prior_diff, error=error),
                      return_posterior_alt(match, mismatchp1, prior_diff=prior_diff, error=error))

    weight = match_prob * kldistance(after_match, before)
    weight += mismatch_prob * kldistance(after_mismatch, before)

    return round(weight, 10)


@lru_cache(maxsize=None)
def return_posterior_ref(match, mismatch, prior_diff=0.001, error=0.01):
    non_error = (1 - error)
    base_is_ref = ((non_error ** match) * (error ** mismatch) * (1 - prior_diff))
    base_is_alt = ((error ** match) * (non_error ** mismatch) * (prior_diff))

    # To prevent division by zero at high coverage
    if (base_is_ref + base_is_alt) == 0:
        ref_alt_prob = sys.float_info.min
    else:
        ref_alt_prob = (base_is_ref + base_is_alt)

    return base_is_ref / ref_alt_prob


@lru_cache(maxsize=None)
def return_posterior_alt(match, mismatch, prior_diff=0.001, error=0.01):
    non_error = (1 - error)
    base_is_ref = ((non_error ** match) * (error ** mismatch) * (1 - prior_diff))
    base_is_alt = ((error ** match) * (non_error ** mismatch) * (prior_diff))

    # To prevent division by zero at high coverage
    if (base_is_ref + base_is_alt) == 0:
        ref_alt_prob = sys.float_info.min
    else:
        ref_alt_prob = (base_is_ref + base_is_alt)

    return base_is_alt / ref_alt_prob


@lru_cache(maxsize=None)
def kldistance(after, before):
    #
    # P || Q (Q to P ) is Sum P(i) (log(Pi) - log(Qi))
    #
    count = 0;
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
    '''
    Helper function that returns True if `i` is an int.
    '''
    try:
        int(i)
        return True
    except ValueError:
        return False


def parse_md_cg_pafline(line):
    """
    This function splits the paf line by \t and returns a named tuple
    """

    PafLine = collections.namedtuple('PafLine', 'read_id chromosome')

    columns = line.split('\t')

    read_id = columns[0]

    ref = columns[5]

    chromosome_temp = ref.split('|')
    chromosome = chromosome_temp[1]

    paffile = PafLine(
        read_id=read_id,
        chromosome=chromosome
    )

    return paffile


def parseMDPAF(pafline):
    t0 = time.time()
    bits = pafline.split()
    sub = "MD:Z:"
    # reflen = int(bits[6])
    # reflen = 20000
    mapstart = int(bits[7])
    mdstr = ""
    for s in filter(lambda x: sub in x, pafline.split()):
        mdstr = s.split(':')[-1]

    '''
    Returns a numpy array with positions for where a `fromMM` base is, in the
    reference genome. This can then be used to intersect with the query string
    to find all desired mutations.
    '''
    # Split MD String at every single [ACGT] or ^:
    mdSub = re.sub(r'([\\^]*[ACGT]+)[0]*', ' \\1 ', mdstr)
    mdSplit = re.split('[ ]+', mdSub)
    # print (mdstr)
    # print (mdSub)
    # print (mdSplit)

    mutArr = np.array([]).astype(int)  # array to hold counts of mismatches
    matArr = np.array([]).astype(int)  # array to hold counts of matches

    # Iterate over Array and replace all mutations from the MD string with the
    # letter of the corresponding reference.
    # eg: 2G1 will produce the numpy array: 'M M G M'
    # All ^[ACGT]* by "D" and the number with a corresponding stretch of "M"
    for i in range(len(mdSplit)):
        mdPos = mdSplit[i]
        if len(mdPos) > 0 and RepresentsInt(mdPos):
            mutArr = np.concatenate((mutArr, np.repeat(0, int(mdPos))))
            matArr = np.concatenate((matArr, np.repeat(1, int(mdPos))))
        elif re.match('\\^', mdPos):
            mutArr = np.concatenate((mutArr, np.repeat(1, len(mdPos) - 1)))
            matArr = np.concatenate((matArr, np.repeat(0, len(mdPos) - 1)))
        elif len(mdPos) == 1:
            mutArr = np.concatenate((mutArr, np.array([1])))
            matArr = np.concatenate((matArr, np.array([0])))
        else:
            # I'm not yet quite sure, if this won't just break at some point.
            # In my BAM/SAM files I have seen rare cases with two consecutive
            # mismatches in the MD tag causing this series of ifs to report incorrect
            # positions if I don't catch this.
            mutArr = np.concatenate((mutArr, np.array([1])))
            matArr = np.concatenate((matArr, np.array([0])))

    # mutArr_length = np.pad(mutArr, (startpad,endpad), 'constant', constant_values = (0,0))
    # matArr_length = np.pad(matArr, (startpad,endpad), 'constant', constant_values = (0,0))

    # Return mismatch positions
    t1 = time.time()
    return mutArr, matArr, mapstart, t1 - t0


def readpaf(lines, paffile='allreadscigar.paf', reflen=5528445):
    paffile = '/Users/mbzros/FAK22471_5ca607ebedbc205cbd837ceb6937de2bfd07e760_0.fastq_output_cigar_md.paf'
    # open file containing example reads
    file = open(paffile, "r")
    # simple counter to limit processing time
    counter = 0
    # line limit
    lines = lines

    # Zero array to hold counts for mismatches and matches of same length as reference
    mismatchholdarray = np.zeros(reflen)
    matchholdarray = np.zeros(reflen)

    # Array to hold information on total coverage of genome
    coveragearray = np.zeros(reflen)

    # Loop through file to read in paf file output and parse MD flag
    for line in file:
        mismatcharray, matcharray, mapstart, runstart = parseMDPAF(line)
        mismatchholdarray[mapstart:mapstart + len(mismatcharray)] += mismatcharray
        matchholdarray[mapstart:mapstart + len(matcharray)] += matcharray
        counter += 1
        # Escape when we have processed max reads for this test.
        if counter >= lines:
            break

    coveragearray = np.sum([matchholdarray, mismatchholdarray], axis=0)

    # print (mismatchholdarray[0:100])
    # print (matchholdarray[0:100])
    # print (coveragearray[0:100])
    return coveragearray, matchholdarray, mismatchholdarray