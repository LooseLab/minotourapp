"""
Functions for calculating the abundance and uncertainty models a la michael
"""
import numpy as np

from .models import CAllValues


def generate_c_all_values(df, task):
    """

    Returns
    -------

    """
    # total_values = (len(targets) - 1) ** 2
    df = df.reset_index().set_index("read_id")
    df["classification_number"] = df["classification_number"].astype(int)
    class_counts_series = df.loc[~df.index.duplicated()].groupby("classification_number")["classification_number"].agg("count")
    # Combine with previous iterations of data
    previous_c_alls = CAllValues.objects.filter(task=task)
    # hold the ids of C-all values that we have saved already so that we don't recreate them below
    saved_class_numbers = set()
    for previous_value in previous_c_alls:
        if previous_value.classification_number in class_counts_series:
            class_counts_series[previous_value.classification_number] += previous_value.classification_count
            # save updated value
            previous_value.classification_count += class_counts_series[previous_value.classification_number]
            previous_value.save()
            saved_class_numbers.add(previous_value.classification_number)
        else:
            class_counts_series.loc[previous_value.classification_number] = previous_value.classification_count
    probabilities_series = class_counts_series/class_counts_series.sum()
    for classification_number, classification_count in class_counts_series.items():
        if classification_number not in saved_class_numbers:
            CAllValues.objects.create(task=task, classification_count=classification_count, classification_number=classification_number)
    return probabilities_series


def generate_c_all(targets, classification_converter, values):
    """
    Generate the c_all values using real data
    Parameters
    ----------
    targets
    classification_converter: dict
    values: pd.core.series.Series
        The values we are inserting into the c_all array from the overarching iterations

    Returns
    -------

    """
    c_all = np.zeros((len(targets), (len(targets) ** 2)))
    c_all_conv = {}
    for i, (k, v) in enumerate(classification_converter.items()):
        c_all_conv[i] = tuple(classification_converter[(tax_id,)] -1 for tax_id in k)
    for k, v in c_all_conv.items():
        value = values.loc[k] if k in values else 0
        np.put(c_all[:, k], v, value)
    # Unclassified is one minus the sum of the previous elements in it's row
    c_all[:, -1] = 1 - c_all[:, :-1].sum(axis=1)
    return c_all


def grammy_fit(classifications, init, draws, species_number, c_all, tolerance=10e-10):
    """
    Grammy model
    Parameters
    ----------
    classifications: numpy.ndarray
        Classifications for each read. Classifications were assigned numbers.
        So if there were two species we were checking for, called them A and B, I set the numbers as follows,
            A only           - 1
            B only           - 2
            A and B        - 3
            Unclassified  - 4
    init: np.ndarray
        init starting point for algorithm - need to explore choice of this
    draws: int
        Number of reads
    species_number: int
        Number of species looking for plus unclassified
    tolerance: float
        Stopping point in the algorithm, its the summed absolute difference between iterations of the algorithm
    c_all: np.ndarray
        the simulated probabilities that determine the probability a read of a certain
        species will receive a certain classification. Each row represents the true species of the read,
        and each column a classification where the column number corresponds to the classification number
        assigned according to the scheme above.

    Returns
    -------
    return start_prob
    """
    err = 1
    start_prob = init
    z = np.zeros((len(classifications), species_number))
    while err > tolerance:
        for i, classification in enumerate(classifications):
            class_prob = c_all[:, classification]
            # todo check with michael what kr, kp, ki represent
            kr = class_prob * init
            zp = kr / np.sum(kr)
            z[i,] = zp
        # @michael todo what do these names mean
        N = z.sum(axis=0)
        pil = N / draws
        err = np.sum(np.absolute(start_prob - pil))
        print(err)
        start_prob = pil
    return start_prob
