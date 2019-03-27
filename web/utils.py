import re

import time

import numpy as np


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

    t0 = time.time()
    bits = line.split()
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
