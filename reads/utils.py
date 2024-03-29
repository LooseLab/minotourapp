import datetime
import itertools

import pandas as pd
from django.db.models import F

from communication.models import Message
from reads.models import (
    RunSummary,
    FastqRead,
    JobMaster,
)


def peek(iterable):
    """
    Peek at the first, return None if the generator is empty
    Parameters
    ----------
    iterable: generator
        Generator to peek at
    Returns
    -------

    """
    try:
        first = next(iterable)
    except StopIteration:
        return None
    return first, itertools.chain([first], iterable)


def create_histogram_series(histogram_value, n50, bin_width, index, name_suffix):
    """
    Return a dictionary of the series values required by the histogram on the live page
    Parameters
    ----------
    histogram_value: int
        The value in bases/estimated bases for that bin
    n50: int
        The n50 value
    bin_width: int
        the width of the bin
    index: int
        The index required
    name_suffix: str
        Suffix to add to the series name. Bases if basecalling enabled else ev


    Returns
    -------

    """
    color = "red" if n50 >= index*bin_width and n50 <= (index+1)*bin_width else "blue"
    name = f"{index*bin_width} - {(index+1)*bin_width} {name_suffix}"
    if color == "red":
        return {"name": name, "color": color, "y": histogram_value, "n50": True}
    return {"name": name, "color": color, "y": histogram_value}


def create_and_save_message(**kwargs):
    """
    Create and save a message to be sent to the user, about a specific action that has happened.
    Parameters
    ----------
    kwargs: dict
        Dictionary of kwargs including, title, user, flowcell, run

    Returns
    -------
    None

    """
    if kwargs["created"]:
        user = kwargs["user"]
        new_message = Message(
            recipient=user,
            sender=user,
            title=kwargs["title"],
            run=kwargs["run"],
            flowcell=kwargs["flowcell"],
        )
        new_message.save()


def return_temp_empty_summary(run):
    """
    If a run has no summary - return an empty RunSummary rather than fail
    :param run: The run that has no summary
    :type run: minknow_data.models.Run
    :return:
    """
    run_summary = RunSummary(
        run=run,
        read_count=0,
        total_read_length=0,
        max_read_length=0,
        avg_read_length=0,
        min_read_length=0,
        first_read_start_time=datetime.datetime.now(datetime.timezone.utc),
        last_read_start_time=datetime.datetime.now(datetime.timezone.utc),
        last_read=0,
        running=False,
    )
    return run_summary


def get_coords(channel, flowcellsize):
    """
    Get the cartesian coordinate on the flowcell layout grid of a channel.
    Parameters
    ----------
    channel: int
        The channel number
    flowcellsize: Total sie of the flowcell

    Returns
    -------

    """
    if flowcellsize == 3000:

        block = (channel - 1) // 250

        remainder = (channel - 1) % 250

        row = remainder // 10

        column = remainder % 10 + block * 10

        return (column, row)

    if flowcellsize == 126:

        return (channel // 8, channel % 8)

    else:

        return flowcell_layout(channel)


def flowcell_layout(channel):
    """
    Lookup dictionary of the channel number to cartesian coordinate of the channel.
    Parameters
    ----------
    channel: int
        The channel number.

    Returns
    -------
    (x, y):tuple
        The x and y coordinate of the channel on the flowcell.

    """

    chanlookup = {
        1: (31, 0),
        2: (31, 1),
        3: (31, 2),
        4: (31, 3),
        5: (31, 4),
        6: (31, 5),
        7: (31, 6),
        8: (31, 7),
        9: (30, 0),
        10: (30, 1),
        11: (30, 2),
        12: (30, 3),
        13: (30, 4),
        14: (30, 5),
        15: (30, 6),
        16: (30, 7),
        17: (29, 0),
        18: (29, 1),
        19: (29, 2),
        20: (29, 3),
        21: (29, 4),
        22: (29, 5),
        23: (29, 6),
        24: (29, 7),
        25: (28, 0),
        26: (28, 1),
        27: (28, 2),
        28: (28, 3),
        29: (28, 4),
        30: (28, 5),
        31: (28, 6),
        32: (28, 7),
        33: (31, 15),
        34: (31, 14),
        35: (31, 13),
        36: (31, 12),
        37: (31, 11),
        38: (31, 10),
        39: (31, 9),
        40: (31, 8),
        41: (30, 15),
        42: (30, 14),
        43: (30, 13),
        44: (30, 12),
        45: (30, 11),
        46: (30, 10),
        47: (30, 9),
        48: (30, 8),
        49: (29, 15),
        50: (29, 14),
        51: (29, 13),
        52: (29, 12),
        53: (29, 11),
        54: (29, 10),
        55: (29, 9),
        56: (29, 8),
        57: (28, 15),
        58: (28, 14),
        59: (28, 13),
        60: (28, 12),
        61: (28, 11),
        62: (28, 10),
        63: (28, 9),
        64: (28, 8),
        65: (3, 0),
        66: (3, 1),
        67: (3, 2),
        68: (3, 3),
        69: (3, 4),
        70: (3, 5),
        71: (3, 6),
        72: (3, 7),
        73: (2, 0),
        74: (2, 1),
        75: (2, 2),
        76: (2, 3),
        77: (2, 4),
        78: (2, 5),
        79: (2, 6),
        80: (2, 7),
        81: (1, 0),
        82: (1, 1),
        83: (1, 2),
        84: (1, 3),
        85: (1, 4),
        86: (1, 5),
        87: (1, 6),
        88: (1, 7),
        89: (0, 0),
        90: (0, 1),
        91: (0, 2),
        92: (0, 3),
        93: (0, 4),
        94: (0, 5),
        95: (0, 6),
        96: (0, 7),
        97: (3, 15),
        98: (3, 14),
        99: (3, 13),
        100: (3, 12),
        101: (3, 11),
        102: (3, 10),
        103: (3, 9),
        104: (3, 8),
        105: (2, 15),
        106: (2, 14),
        107: (2, 13),
        108: (2, 12),
        109: (2, 11),
        110: (2, 10),
        111: (2, 9),
        112: (2, 8),
        113: (1, 15),
        114: (1, 14),
        115: (1, 13),
        116: (1, 12),
        117: (1, 11),
        118: (1, 10),
        119: (1, 9),
        120: (1, 8),
        121: (0, 15),
        122: (0, 14),
        123: (0, 13),
        124: (0, 12),
        125: (0, 11),
        126: (0, 10),
        127: (0, 9),
        128: (0, 8),
        129: (7, 0),
        130: (7, 1),
        131: (7, 2),
        132: (7, 3),
        133: (7, 4),
        134: (7, 5),
        135: (7, 6),
        136: (7, 7),
        137: (6, 0),
        138: (6, 1),
        139: (6, 2),
        140: (6, 3),
        141: (6, 4),
        142: (6, 5),
        143: (6, 6),
        144: (6, 7),
        145: (5, 0),
        146: (5, 1),
        147: (5, 2),
        148: (5, 3),
        149: (5, 4),
        150: (5, 5),
        151: (5, 6),
        152: (5, 7),
        153: (4, 0),
        154: (4, 1),
        155: (4, 2),
        156: (4, 3),
        157: (4, 4),
        158: (4, 5),
        159: (4, 6),
        160: (4, 7),
        161: (7, 15),
        162: (7, 14),
        163: (7, 13),
        164: (7, 12),
        165: (7, 11),
        166: (7, 10),
        167: (7, 9),
        168: (7, 8),
        169: (6, 15),
        170: (6, 14),
        171: (6, 13),
        172: (6, 12),
        173: (6, 11),
        174: (6, 10),
        175: (6, 9),
        176: (6, 8),
        177: (5, 15),
        178: (5, 14),
        179: (5, 13),
        180: (5, 12),
        181: (5, 11),
        182: (5, 10),
        183: (5, 9),
        184: (5, 8),
        185: (4, 15),
        186: (4, 14),
        187: (4, 13),
        188: (4, 12),
        189: (4, 11),
        190: (4, 10),
        191: (4, 9),
        192: (4, 8),
        193: (11, 0),
        194: (11, 1),
        195: (11, 2),
        196: (11, 3),
        197: (11, 4),
        198: (11, 5),
        199: (11, 6),
        200: (11, 7),
        201: (10, 0),
        202: (10, 1),
        203: (10, 2),
        204: (10, 3),
        205: (10, 4),
        206: (10, 5),
        207: (10, 6),
        208: (10, 7),
        209: (9, 0),
        210: (9, 1),
        211: (9, 2),
        212: (9, 3),
        213: (9, 4),
        214: (9, 5),
        215: (9, 6),
        216: (9, 7),
        217: (8, 0),
        218: (8, 1),
        219: (8, 2),
        220: (8, 3),
        221: (8, 4),
        222: (8, 5),
        223: (8, 6),
        224: (8, 7),
        225: (11, 15),
        226: (11, 14),
        227: (11, 13),
        228: (11, 12),
        229: (11, 11),
        230: (11, 10),
        231: (11, 9),
        232: (11, 8),
        233: (10, 15),
        234: (10, 14),
        235: (10, 13),
        236: (10, 12),
        237: (10, 11),
        238: (10, 10),
        239: (10, 9),
        240: (10, 8),
        241: (9, 15),
        242: (9, 14),
        243: (9, 13),
        244: (9, 12),
        245: (9, 11),
        246: (9, 10),
        247: (9, 9),
        248: (9, 8),
        249: (8, 15),
        250: (8, 14),
        251: (8, 13),
        252: (8, 12),
        253: (8, 11),
        254: (8, 10),
        255: (8, 9),
        256: (8, 8),
        257: (15, 0),
        258: (15, 1),
        259: (15, 2),
        260: (15, 3),
        261: (15, 4),
        262: (15, 5),
        263: (15, 6),
        264: (15, 7),
        265: (14, 0),
        266: (14, 1),
        267: (14, 2),
        268: (14, 3),
        269: (14, 4),
        270: (14, 5),
        271: (14, 6),
        272: (14, 7),
        273: (13, 0),
        274: (13, 1),
        275: (13, 2),
        276: (13, 3),
        277: (13, 4),
        278: (13, 5),
        279: (13, 6),
        280: (13, 7),
        281: (12, 0),
        282: (12, 1),
        283: (12, 2),
        284: (12, 3),
        285: (12, 4),
        286: (12, 5),
        287: (12, 6),
        288: (12, 7),
        289: (15, 15),
        290: (15, 14),
        291: (15, 13),
        292: (15, 12),
        293: (15, 11),
        294: (15, 10),
        295: (15, 9),
        296: (15, 8),
        297: (14, 15),
        298: (14, 14),
        299: (14, 13),
        300: (14, 12),
        301: (14, 11),
        302: (14, 10),
        303: (14, 9),
        304: (14, 8),
        305: (13, 15),
        306: (13, 14),
        307: (13, 13),
        308: (13, 12),
        309: (13, 11),
        310: (13, 10),
        311: (13, 9),
        312: (13, 8),
        313: (12, 15),
        314: (12, 14),
        315: (12, 13),
        316: (12, 12),
        317: (12, 11),
        318: (12, 10),
        319: (12, 9),
        320: (12, 8),
        321: (19, 0),
        322: (19, 1),
        323: (19, 2),
        324: (19, 3),
        325: (19, 4),
        326: (19, 5),
        327: (19, 6),
        328: (19, 7),
        329: (18, 0),
        330: (18, 1),
        331: (18, 2),
        332: (18, 3),
        333: (18, 4),
        334: (18, 5),
        335: (18, 6),
        336: (18, 7),
        337: (17, 0),
        338: (17, 1),
        339: (17, 2),
        340: (17, 3),
        341: (17, 4),
        342: (17, 5),
        343: (17, 6),
        344: (17, 7),
        345: (16, 0),
        346: (16, 1),
        347: (16, 2),
        348: (16, 3),
        349: (16, 4),
        350: (16, 5),
        351: (16, 6),
        352: (16, 7),
        353: (19, 15),
        354: (19, 14),
        355: (19, 13),
        356: (19, 12),
        357: (19, 11),
        358: (19, 10),
        359: (19, 9),
        360: (19, 8),
        361: (18, 15),
        362: (18, 14),
        363: (18, 13),
        364: (18, 12),
        365: (18, 11),
        366: (18, 10),
        367: (18, 9),
        368: (18, 8),
        369: (17, 15),
        370: (17, 14),
        371: (17, 13),
        372: (17, 12),
        373: (17, 11),
        374: (17, 10),
        375: (17, 9),
        376: (17, 8),
        377: (16, 15),
        378: (16, 14),
        379: (16, 13),
        380: (16, 12),
        381: (16, 11),
        382: (16, 10),
        383: (16, 9),
        384: (16, 8),
        385: (23, 0),
        386: (23, 1),
        387: (23, 2),
        388: (23, 3),
        389: (23, 4),
        390: (23, 5),
        391: (23, 6),
        392: (23, 7),
        393: (22, 0),
        394: (22, 1),
        395: (22, 2),
        396: (22, 3),
        397: (22, 4),
        398: (22, 5),
        399: (22, 6),
        400: (22, 7),
        401: (21, 0),
        402: (21, 1),
        403: (21, 2),
        404: (21, 3),
        405: (21, 4),
        406: (21, 5),
        407: (21, 6),
        408: (21, 7),
        409: (20, 0),
        410: (20, 1),
        411: (20, 2),
        412: (20, 3),
        413: (20, 4),
        414: (20, 5),
        415: (20, 6),
        416: (20, 7),
        417: (23, 15),
        418: (23, 14),
        419: (23, 13),
        420: (23, 12),
        421: (23, 11),
        422: (23, 10),
        423: (23, 9),
        424: (23, 8),
        425: (22, 15),
        426: (22, 14),
        427: (22, 13),
        428: (22, 12),
        429: (22, 11),
        430: (22, 10),
        431: (22, 9),
        432: (22, 8),
        433: (21, 15),
        434: (21, 14),
        435: (21, 13),
        436: (21, 12),
        437: (21, 11),
        438: (21, 10),
        439: (21, 9),
        440: (21, 8),
        441: (20, 15),
        442: (20, 14),
        443: (20, 13),
        444: (20, 12),
        445: (20, 11),
        446: (20, 10),
        447: (20, 9),
        448: (20, 8),
        449: (27, 0),
        450: (27, 1),
        451: (27, 2),
        452: (27, 3),
        453: (27, 4),
        454: (27, 5),
        455: (27, 6),
        456: (27, 7),
        457: (26, 0),
        458: (26, 1),
        459: (26, 2),
        460: (26, 3),
        461: (26, 4),
        462: (26, 5),
        463: (26, 6),
        464: (26, 7),
        465: (25, 0),
        466: (25, 1),
        467: (25, 2),
        468: (25, 3),
        469: (25, 4),
        470: (25, 5),
        471: (25, 6),
        472: (25, 7),
        473: (24, 0),
        474: (24, 1),
        475: (24, 2),
        476: (24, 3),
        477: (24, 4),
        478: (24, 5),
        479: (24, 6),
        480: (24, 7),
        481: (27, 15),
        482: (27, 14),
        483: (27, 13),
        484: (27, 12),
        485: (27, 11),
        486: (27, 10),
        487: (27, 9),
        488: (27, 8),
        489: (26, 15),
        490: (26, 14),
        491: (26, 13),
        492: (26, 12),
        493: (26, 11),
        494: (26, 10),
        495: (26, 9),
        496: (26, 8),
        497: (25, 15),
        498: (25, 14),
        499: (25, 13),
        500: (25, 12),
        501: (25, 11),
        502: (25, 10),
        503: (25, 9),
        504: (25, 8),
        505: (24, 15),
        506: (24, 14),
        507: (24, 13),
        508: (24, 12),
        509: (24, 11),
        510: (24, 10),
        511: (24, 9),
        512: (24, 8),
    }

    return chanlookup[channel]


def getn50(lens):
    """
    Get the n50 of a read set
    Parameters
    ----------
    lens: list
        List of read lengths

    Returns
    -------
    l: int
        The length of the read that is the n50 length
    """
    # TODO could be rewritten into numpy quite easily
    h = sum(lens) / 2
    t = 0
    for l in lens:
        t += l
        if t >= h:
            return l


def pause_job(job_master):
    """
    Pause a job master.
    Parameters
    ----------
    job_master: reads.models.JobMaster
        The JobMaster to pause.
    Returns
    -------
    job_master, return_message: (reads.models.JobMaster, str)
        The modified job master object and the return message
    """

    if job_master.paused:
        job_master.paused = False
        paused_status = "un-paused"
    else:
        job_master.paused = True
        paused_status = "paused"

    job_master.save()
    return_message = f"Successfully {paused_status} {job_master.job_type.name} task, id: {job_master.id}"
    return job_master, return_message


def clear_artic_command_job_masters(flowcell_id):
    """
    Clear the non outdated JobMasters that ran commands
    Parameters
    ----------
    flowcell_id: int
        The primary key of the flowcell that has had the task reset

    Returns
    -------
    None
    """
    JobMaster.objects.filter(flowcell_id=flowcell_id, job_type_id=17).delete()


def get_fastq_df(desired_yield, flowcell_pk, avg_read_len, task):
    """
    Get a set yield of fastq from the database, and return as a dataframe

    Parameters
    ----------
    desired_yield: int
        Desired yield in megabases
    flowcell_pk: int
        The flowcell primary key
    avg_read_len: int
        The average read length in bases
    task: reads.models.JobMaster
        The task that we are fetching the fastq for

    Returns
    -------

    """
    desired_yield = desired_yield * 1000000
    chunk_size = round(desired_yield / avg_read_len)
    print(f"Fetching reads in chunks of {chunk_size} for alignment.")
    fasta_df_barcode = pd.DataFrame().from_records(
        FastqRead.objects.filter(
            flowcell_id=flowcell_pk, id__gt=task.last_read
        )[:chunk_size].values(
            "read_id",
            "barcode_name",
            "sequence",
            "id",
            "run_id",
            "type__name",
            "run_id",
            "barcode_id",
            "is_pass",
            read_type_id=F("type_id"),
        )
    )
    if not fasta_df_barcode.empty:
        last_read = fasta_df_barcode.tail(1).iloc[0].id
    if fasta_df_barcode.empty:
        print("No fastq found. Skipping this iteration.")
        task.running = False
        task.save()
        return None, None, pd.DataFrame()
    read_count = fasta_df_barcode.shape[0]
    return read_count, last_read, fasta_df_barcode
