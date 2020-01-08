import numpy as np
import pandas as pd
from pathlib import Path
from collections import namedtuple
from django.conf import settings


def make_results_directory(flowcell_id, task_id):
    """
    Make directory to store results of this task
    :param flowcell_id: The ID of the flowcell that produced these fastqs
    :param task_id: The ID of the Jobmaster that corresponds to this EB task
    :return: The path into this directory
    """
    current_working_directory = Path.cwd()
    # current_working_directory = getattr(settings, "REFERENCE_LOCATION", None)

    results_dir = Path(
        f"{current_working_directory}/readuntil/Temp_results/EB{flowcell_id}_{task_id}"
    )


    if not results_dir.exists():
        Path.mkdir(results_dir)
    return results_dir


def setup_chromosome(reference_length):
    return np.zeros(reference_length)


def add_chromosome(s_store, reference, reference_length):
    s_store[reference] = dict()
    for thing in ["I", "IC", "D", "A", "T", "G", "C", "M"]:
        s_store[reference][thing] = setup_chromosome(int(reference_length))
    return s_store


def sum_dictionaries(s_store, d, reference, mapstart, mapend):
    # dictkeys = d.keys()
    for dictkey in d.keys():

        s_store[reference][dictkey][mapstart:mapend] += d[dictkey]
    return s_store


def _rstrip_tuple(s, chars):
    """
    Return a tuple of the initial string and the chars removed from the right
    """
    h = s.rstrip(chars)
    return h, s[len(h) :]


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
        yield _conv_type(e, int)


def parse_MD(s, match=0, mismatch=1):
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


def _expand_dict_in_series(df, field):
    """Convert a Series of dict to Series and add to the original DataFrame

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with a Series of dict
    field : str
        The Series of dicts to expand

    Returns
    -------
    pd.DataFrame
        The orignal DataFrame with extra Series from the dicts
    """
    return df.join(pd.DataFrame(df.pop(field).tolist()))


def _parse_tags(tags):
    """Convert a list of SAM style tags, from a PAF file, to a dict

    Parameters
    ----------
    tags : list
        A list of SAM style tags

    Returns
    -------
    dict
        Returns dict of SAM style tags
    """
    c = {"i": int, "A": str, "f": float, "Z": str}
    return {
        key: _conv_type(val, c[tag])
        for key, tag, val in (x.split(":") for x in tags)
    }


def _conv_type(s, func):
    """Generic converter, to change strings to other types

    Parameters
    ----------
    s : str
        A string that represents another type
    func : type
        type to apply to s

    Returns
    -------
    The type of func, otherwise str
    """
    try:
        return func(s)
    except ValueError:
        return s


def _format_records(record):
    """Helper function to make fields the right type
    """
    return [_conv_type(x, int) for x in record]


def _paf_generator(file_like, fields=None):
    """Generator that returns namedtuples from a PAF file

    Parameters
    ----------
    file_like : file-like object
        File-like object
    fields : list
        List of field names to use for records, must have 13 entries.

    Yields
    ------
    namedtuple
        Correctly formatted PAF record and a dict of extra tags

    Raises
    ------
    ValueError
    """
    if len(fields) != 13:
        raise ValueError("{} fields provided, expected 13".format(len(fields)))
    PAF = namedtuple("PAF", fields)
    for record in file_like:
        record = record.strip().split("\t")
        yield PAF(*_format_records(record[:12]), _parse_tags(record[12:]))


def parse_PAF(file_like, fields=None, dataframe=False):
    """Read a minimap2 PAF file as either an iterator or a pandas.DataFrame
    Parameters
    ----------
    file_like : file-like object
        Object with a read() method, such as a file handler or io.StringIO.
    fields : list
        List of field names to use for records, must have 13 entries. Default:
        ["read_name", "query_length", "query_start", "query_end", "strand", "target_name", "target_length",
        "target_start", "target_end", "residue_matches", "alignment_block_length", "mapping_quality", "tags"]
    dataframe : bool
        When True a pandas.DataFrame is returned with Series named as the `fields` parameter.
        SAM tags are expanded into Series as well and given their specified types.
    Returns
    -------
    iterator or pandas.DataFrame when dataframe is True
    """
    if fields is None:
        fields = [
            "read_name",
            "query_length",
            "query_start",
            "query_end",
            "strand",
            "target_name",
            "target_length",
            "target_start",
            "target_end",
            "residue_matches",
            "alignment_block_length",
            "mapping_quality",
            "tags",
        ]
    if dataframe:
        return _expand_dict_in_series(
            pd.DataFrame(_paf_generator(file_like, fields=fields)), fields[-1]
        )
    else:
        return _paf_generator(file_like, fields=fields)


def readfq(fp):
    """
    DocString
    """
    # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        desc, name, seqs, last = last[1:], last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield desc, name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield desc, name, seq, "".join(
                        seqs
                    )  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield desc, name, seq, None  # yield a fasta record instead
                break


def _lstrip_tuple(s, chars):
    """Return a tuple of the initial string and the chars removed from the left
    """
    remains = s.lstrip(chars)
    return remains, s[: len(s) - len(remains)]


def cg_generator(s):
    """Given a CIGAR string yields elements in order

    Parameters
    ----------
    s : str
        A "cg:Z:" string from a mapping without the tag

    Returns
    -------
    tuple
        Generator that partions the CIGAR tag in order,
        yielding a tuple of the count and the operation
    """
    while s:
        s, c = _lstrip_tuple(s, "0123456789")
        s, o = _lstrip_tuple(s, "MIDNSHPX=")
        yield _conv_type(c, int), _conv_type(o, int)


def parse_cigar(cigar_string, read_string, chunk_length):
    """Return a dictionary of arrays

    Parameters
    ----------
    cigar_string : str
        A CIGAR string
    read_string : str
        Read that the CIGAR string describes
    chunk_length : int
        Length of the reference covered by the alignment.

    Returns
    -------
    array
    dict
    """
    d = {}
    # Array of insertions
    d["I"] = np.zeros(chunk_length)
    # Array of insertion counts
    d["IC"] = np.zeros(chunk_length)
    query_pos = 0
    query_array = []
    read_bases = tuple(read_string)

    for count, operation in cg_generator(cigar_string):
        # print (operation,count)
        # Not aligned read section
        if operation is "S":
            # print ("skipping")
            pass

        # Not a deletion or insertion. Its 0:M
        elif operation is "M":
            # Extend query_array with bases from read
            # print ("match")
            query_array.extend(read_bases[query_pos: query_pos + count])

        elif operation is "I":
            # print ("insertion")
            # Extend query_array with bases from read
            # query_array.extend(read_bases[query_pos:query_pos + count])

            d["I"][query_pos] += 1
            d["IC"][query_pos] += count
            count = 1

        elif operation is "D":
            # print ("deletion")
            query_array.extend(["D"] * count)
            continue
        else:
            print("Operation: '{}' not accounted for.".format(operation))
        # Increment query position by count
        # TODO IS THIS BELOW CODE COOL
        if operation is not "I":
            query_pos += count

    query_array = np.fromiter(query_array, "U1")
    for base in np.unique(query_array):
        d[base] = np.zeros(len(query_array))
        np.add.at(d[base], np.argwhere(query_array == base), 1)
    return d


def multi_array_results(ref_length):
    """
    Create a structured array to contain the counts of the nucleotide bases mapped to each position on the reference
    :param ref_length: The length of the reference chromosome
    :type ref_length: int
    :return: A filled multi dimensional array, the length of the reference chromosome
    """
    # multiply out a list to create as many copies as there are bases
    a = [tuple(0 for i in range(9))] * ref_length
    # Create a structure data array with named data types - create 9 copies of the above list in a list,
    # and make the structured array with them
    multi_arr = np.array(
        [a for i in range(8)][0],
        dtype=[
            ("A", np.uint16),
            ("C", np.uint16),
            ("G", np.uint16),
            ("T", np.uint16),
            ("D", np.uint16),
            ("I", np.uint16),
            ("IC", np.uint16),
            ("M", np.uint16),
            ("U", np.bool),
        ],
    )
    # return our new numpy array
    return multi_arr


def initialise_priors(
    reference, genotypes, theta=0.009, r=0.11
):
    """

    :param reference:
    :param genotypes:
    :param theta:
    :param r:
    :return:
    """
    if theta * (1.0 + r) > 1.0:
        print(
            "PROBLEM: theta too high or indel rate too high, almost all"
            " genome expected to be different from the reference"
        )
        return
    # initialise the empty list that will hold our resulting priors
    priors = []
    # Create an array where each letter of the reference is an element in the array
    query_array = np.fromiter(reference, "U1")
    # each element in all bu the last genotypes element
    for g in genotypes[: len(genotypes) - 1]:
        # append an array to priors for wach of the four base types, where if
        #  the reference base is the same as the for loop base
        # the element is the same as 1 - theta * 1 = r, if it's not the same,
        #  the element is the same as theta divided by 5-2
        priors.append(
            np.where(
                (query_array == g),
                (1.0 - (theta * (1.0 + r))),
                theta / (len(genotypes) - 2),
            )
        )
    # add one final array to the priors list where there is as many elements as
    #  reference bases and the element value is theta * r
    priors.append(np.full(query_array.size, theta * r))
    # turn the priors list into an numpy array
    priors = np.array(priors)
    # return the priors array transformed so it's as many rows as the reference
    # with 5 element in each, one from each array above
    return priors.T


def initialise_posteriors_mod(priors):
    # just copy the array, as they are the same
    posteriors = priors.copy()
    # return it
    return posteriors


def position_benefit(posterior_array, phi):
    # Shannon entropy of the current posterior distribution
    entropies = -np.sum(
        posterior_array * np.log(posterior_array, where=posterior_array != 0),
        axis=1,
    )
    # probability of observing a new base of each type at the considered position
    observation_probabilities = np.array([])
    # A list to append arrays of new posteriors to
    new_posteriors = list()
    # For each column in the 2d array
    for i in range(len(phi)):
        # Add the sum of each phi column * the sum of the row of the 2d posterior array
        observation_probabilities = np.concatenate(
            (
                observation_probabilities,
                np.sum(phi[:, i] * posterior_array, axis=1),
            )
        )
        # New posteriors is the same as the observation probabilities, but each element in a row is normalised
        #  by the sum of the row
        # So E is the posterior array mulitplied by a column of phi
        e = phi[:, i] * posterior_array

        e_sum = e.sum(axis=1, keepdims=True)
        # Now divide each row of e by it's sum, should be one for each posterior, being careful not to divide by 0
        e = np.divide(e, e_sum, where=e_sum != 0)
        # Append the normalise rows to this list
        new_posteriors.append(e)

    # Reshape the probabilitiesto the len of phi columns and number of posteriors rows and transfrom it
    observation_probabilities = np.resize(
        observation_probabilities, (len(phi), len(posterior_array))
    ).T

    new_entropies = list()
    # For each phi, so 5
    for i in range(len(phi)):
        # Create a new Entropies list that contains the negative sum value of new posteriors * observation_probabilites,
        # ravelled to 1d, then shaped into a 2d array (so the multiplication of the two arrays is respective),
        # multiplied by the log of the new log of the posteriors
        new_entropies.append(
            -np.sum(
                new_posteriors[i]
                * observation_probabilities[:, i, None].ravel()[:, np.newaxis]
                * np.log(new_posteriors[i], where=new_posteriors[i] != 0),
                axis=1,
            )
        )

    new_entropies = np.array(new_entropies)

    new_entropies[new_entropies == -0] = 0

    # Sum the columns of the new_entropies array, returns 1d array
    new_entropies = new_entropies.sum(axis=0)

    # Return the difference between the two entropies
    return entropies - new_entropies


def update_posteriors(
    test_posteriors, test_counts, testphi
):
    # Check your phis aren't too high
    # test_posteriors array is a nested array the length of the reference,
    # there are 5 elements in each nest, one for each base plus deletions
    # test_counts is the counts structured data array, as many elements
    # as there are reference bases, 9 elements in each array
    # test phi, is the phis

    # Get total count of bases
    reference_length = test_posteriors.shape[0]

    # Expand test counts in 2d array, with 5 copies of each of each test count element
    d = np.resize(
        np.repeat(testphi[None, ...], 5, axis=0),
        (reference_length * 5, 5),
    )
    # Expand test counts in 2d array, with 5 copies of each of each test count element

    # Do in place exponentials, memory efficient
    d **= np.repeat(
        np.array(test_counts.tolist()),
        test_posteriors.shape[1],
        axis=0,
    )

    # Get all the posteriors as 1d array
    # insert_test_posteriors = np.resize(np.repeat(test_posteriors, 5), (reference_length*5,5)).T[0]
    insert_test_posteriors = test_posteriors.ravel()

    # Create an array of zeros with as many rows as there positions, with 6 elements in each row
    arr = np.zeros((d.shape[0], 6))
    # Set all but the first row as the exponential multiplication of the counts and phis
    arr[:, 1:] = d
    # Set the first column as the test posteriors, one value to each row
    arr[:, 0] = insert_test_posteriors

    results = np.cumprod(arr, axis=1)
    # Ravel the last column of the 2d results array
    results = results[:, -1:].ravel()

    # Np resize the results to be a 2d array of 5 elements, with as many rows as you have reference positions,
    #  so correct shape as beginning posteriors
    results = np.resize(results, (reference_length, 5))
    # print(results)
    # Normalise the posteriors against the sum of the row, taking care to return 0 if we have a divide by 0
    posteriors = np.divide(
        results,
        results.sum(axis=1)[:, None],
        out=np.zeros_like(results),
        where=results.sum(axis=1)[:, None] != 0,
    )

    # posteriors.dtype = np.float16

    return posteriors


def calculate_read_benefits_fixed_read_length(benefits, lam, m):
    # Forward Benefits
    # Number of reference positions, as there is one element in benefits array for each position
    N = benefits.shape[0]
    # Get the partials value to start with, sum of the first 4 benfits values
    partialS = benefits[: min(lam, N)].sum()
    # Get the number of elements that would have been in Nics for loop
    benefits_test = benefits[: N - m]
    # Add the first partials into start of array, so we can use it for cumulative sum
    benefits_test = np.insert(
        benefits_test, 0, partialS, axis=0
    )
    # Make all but the first element it's negative, as it is cumulative subtraction
    benefits_test[1:] = -benefits_test[1:]
    # Cumulative sum the array
    partials_1 = np.cumsum(benefits_test)
    # Return array to original negative/positive
    benefits_test[1:] = -benefits_test[1:]
    # Get the values that would have met the if condition in nics function
    cond_bens = benefits[N - (N - lam):]
    print(f"N is {N}")
    print(f"m is {m}")
    print(f"benefits_test size is {benefits_test.size}")
    print(f"cond bens size is {cond_bens.size}")
    print(f"partilas 1 size is {partials_1.size}")
    # Add the cumulative sum array of those values where they would have been added in for loop
    partials_1[1 : 1 + cond_bens.size] += np.cumsum(cond_bens)
    # Add the sum of those values to every element after where they would have been to make up for their absence
    partials_1[1 + cond_bens.size :] += cond_bens.sum()
    # rename to results array
    read_benefits_forward = partials_1
    # Get the hard truncation end value of 0
    zero_benefit = np.zeros(m - 1, dtype=np.int8)
    # Add the 0 array onto the end of the results array!
    read_benefits_forward = np.concatenate(
        (read_benefits_forward, zero_benefit), axis=0
    )

    # ####### READ BENEFITS REVERSE ########
    # Create zeros for hard cut off at start, not enough read length?
    read_benefits_rev = np.zeros(m - 1)
    # sum the elements that proceed m
    partialS_rev = benefits[: m - 1].sum()
    # Create the array that contains the elemnts that would have been in the for loop
    bens_rev = benefits[m - 1:]
    # Prepend our partialS value to the start of the array
    bens_rev = np.insert(
        bens_rev, 0, partialS_rev, axis=0
    )
    # Cumulatively sum the array
    partials1_rev = np.cumsum(bens_rev)
    # Drop the first value
    partials1_rev = partials1_rev[1:]
    # Get the values that met the condition in Nics 2nd for loop, excluding the m-1 zeros we have at the front
    cond_bens_rev = benefits[lam - m + 1 + (m - 1):]
    # Cumulataively sum them
    cond_bens_rev = np.cumsum(cond_bens_rev)
    # Subtract them from the array at the correct postion (In this case at the end)
    partials1_rev[
        lam - m + 1 : lam - m + 1 + cond_bens_rev.size
    ] -= cond_bens_rev
    # Concat our results to the results array
    read_benefits_rev = np.concatenate(
        (read_benefits_rev, partials1_rev), axis=0
    )

    return read_benefits_forward, read_benefits_rev


def calculate_costs_fixed_read_length(N, rho, m, lam):
    # Create empty array of zeros correct length
    read_costs_forward = np.zeros(N)
    # Add lam to it where lam is smaller than N - index
    read_costs_forward[:lam] += lam
    # The value of N minus the index where it is smaller than lam
    read_costs_forward[-lam + 1 :] += np.arange(lam - 1, 0, -1)
    # Where m is smaller than the reference minus the element index, subtract m+rho
    read_costs_forward[:-m] -= m + rho
    # Where m is larger than the reference minus the element index, subtract the reference lenght - index +rho
    read_costs_forward[-m:] -= np.arange(m, 0, -1) + rho
    # Reverse it for the costs
    read_costs_reverse = read_costs_forward[::-1].copy()

    return read_costs_forward, read_costs_reverse


def calculate_scores(
    read_benefits_forward,
    read_benefits_reverse,
    m_benefits_forward,
    m_benefits_reverse,
    read_costs_forward,
    read_costs_reverse,
):
    # Create an array of correct length, with the length of the reference
    N = read_benefits_forward.size
    scores_F = np.zeros(N)
    # Just list infinity if the cost is negative or too small
    index = read_costs_forward <= 0.0000001
    scores_F[index] = np.inf
    # If it is bigger than 0.0000001 do that simple calculation

    scores_F[~index] = (
                               read_benefits_forward[~index] - m_benefits_forward[~index]
    ) / read_costs_forward[~index]
    # Do the same for reverse
    scores_R = np.zeros(N)
    # Infinity for something that's too small
    index = read_costs_reverse <= 0.0000001
    scores_R[index] = np.inf
    # Do the same calculation as above with the reverse values
    scores_R[~index] = (
                               read_benefits_reverse[~index] - m_benefits_reverse[~index]
    ) / read_costs_reverse[~index]

    return scores_F, scores_R


def find_strategy_uniform(
    scores_f,
    scoresR,
    readBenefitsF,
    readBenefitsR,
    mBenefitsF,
    mBenefitsR,
    costsF,
    costsR,
    rho,
    m,
    alpha,
):
    acceptF = []
    acceptR = []
    totScores = []
    # 1 -  sort the scores and initialise the strategy to always reject
    for pos in range(len(scores_f)):
        totScores.append((scores_f[pos], 1, pos))
        acceptF.append(0)
    for pos in range(len(scoresR)):
        totScores.append((scoresR[pos], 0, pos))
        acceptR.append(0)
    totScores.sort(reverse=True)
    # 2 - initialize the strategy benefit and cost
    StrategyBenefit = sum(mBenefitsF) + sum(mBenefitsR)
    StrategyCost = (alpha + m + rho) * (len(scores_f) + len(scoresR)) - m * (
        m - 1
    )  # the last bit is to account for positions very close to the end, for which lamda_i<m .
    #  Maybe I should assume 0 benefit for these?
    # 3 - find the best strategy
    for score in totScores:
        if score[1] == 1:  # if it's a forward read
            newB = (
                StrategyBenefit
                + readBenefitsF[score[2]]
                - mBenefitsF[score[2]]
            )
            newC = StrategyCost + costsF[score[2]]
            if newB / newC >= StrategyBenefit / StrategyCost:
                acceptF[score[2]] = 1
                StrategyBenefit = newB
                StrategyCost = newC
            else:
                break
        else:  # reverse read
            newB = (
                StrategyBenefit
                + readBenefitsR[score[2]]
                - mBenefitsR[score[2]]
            )
            newC = StrategyCost + costsR[score[2]]
            if newB / newC >= StrategyBenefit / StrategyCost:
                acceptR[score[2]] = 1
                StrategyBenefit = newB
                StrategyCost = newC
            else:
                break

    return np.array(acceptF, dtype=np.bool), np.array(acceptR, dtype=np.bool)
