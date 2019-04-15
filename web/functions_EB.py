# ##Function modified from https://raw.githubusercontent.com/lh3/readfq/master/readfq.py

import matplotlib.pyplot as plt
import numpy as np
import math
from functools import lru_cache
import re, sys, os
import time
import gzip
import pandas as pd

from collections import namedtuple

from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)


# Set up the logger to write to logging file


def setup_chromosome(reference_length):
    return np.zeros(reference_length)


def add_chromosome(s_store, reference, reference_length):
    s_store[reference] = dict()
    for thing in ['I', 'IC', 'D', "A", "T", "G", "C", "M"]:
        s_store[reference][thing] = setup_chromosome(int(reference_length))
    return s_store


def sum_dictionaries(s_store, d, reference, mapstart, mapend):
    # dictkeys = d.keys()
    for dictkey in d.keys():
        # print (dictkey,len(d[dictkey]))
        s_store[reference][dictkey][mapstart:mapend] += d[dictkey]
    return s_store


def _rstrip_tuple(s, chars):
    """
    Return a tuple of the initial string and the chars removed from the right
    """
    h = s.rstrip(chars)
    return h, s[len(h):]


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
        raise ValueError(
            "{} fields provided, expected 13".format(
                len(fields)
            )
        )
    PAF = namedtuple("PAF", fields)
    for record in file_like:
        logger.info(f"paf line is {record}")
        record = record.strip().split("\t")
        yield PAF(
            *_format_records(record[:12]),
            _parse_tags(record[12:])
        )


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
        When True a pandas.DataFrame is returned with Series named as the `fields` parameter. SAM tags are expanded
        into Series as well and given their specified types.

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
            pd.DataFrame(
                _paf_generator(file_like, fields=fields)
            ),
            fields[-1],
        )
    else:
        return _paf_generator(file_like, fields=fields)


def _lstrip_tuple(s, chars):
    """Return a tuple of the initial string and the chars removed from the left
    """
    remains = s.lstrip(chars)
    return remains, s[:len(s) - len(remains)]


def _cg_generator(s):
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


def parse_CIGAR(cigar_string, read_string, alignment_length):
    """Return an array and dictionary

    Parameters
    ----------
    cigar_string : str
        A CIGAR string
    read_string : str
        Read that the CIGAR string describes
    alignment_length : int
        Length of the alignment

    Returns
    -------
    array
    dict
    """
    # Array of insertions
    # Array of insertion counts
    d = {"I": np.zeros(alignment_length),
         "IC": np.zeros(alignment_length)}

    query_pos = 0

    query_array = []

    read_bases = tuple(read_string)

    for count, operation in _cg_generator(cigar_string):
        # Not aligned read section
        if operation is 'S':
            pass

        # Not a deletion or insertion. Its 0:M
        elif operation is 'M':
            # Extend query_array with bases from read
            query_array.extend(read_bases[query_pos:query_pos + count])

        elif operation is 'I':
            # Extend query_array with bases from read
            # query_array.extend(read_bases[query_pos:query_pos + count])
            query_array_len = len(query_array) - 1
            d['I'][query_array_len] += 1
            d["IC"][query_array_len] += count
        #             count=1
        elif operation is 'D':
            query_array.extend(["D"] * count)
        else:
            print("Operation: '{}' not accounted for.".format(operation))
        # Increment query position by count
        if operation is not "D":
            query_pos += count

    query_array = np.fromiter(query_array, "U1")
    for base in np.unique(query_array):
        d[base] = np.zeros(len(query_array))
        np.add.at(d[base], np.argwhere(query_array == base), 1)
    return query_array, d
