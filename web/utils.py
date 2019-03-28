import re

import time

import numpy as np
import collections


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