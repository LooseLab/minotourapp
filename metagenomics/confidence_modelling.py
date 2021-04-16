"""
Functions for calculating the abundance and uncertainty models a la michael
"""
import numpy as np
from scipy.linalg import lapack

from .models import CAllValues, EstimatedAbundance, UncertaintyProbability

class_probs = np.array(
    [
        [8.3003e-01, 2.0000e-05],
        [0.0000e00, 6.3356e-01],
        [0.0000e00, 0.0000e00],
        [1.6997e-01, 3.6642e-01],
    ]
).T


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


def save_abundances(tax_id_2_prob, names, task, class_conv):
    """
    Save the abundances by tax_id
    Parameters
    ----------
    tax_id_2_prob: dict
        Dictionary keyed tax_id to estimated abundance value
    names: dict
        tax_ids keyed to names
    task: reads.models.JobMaster
        The task for this metagenomics task
    class_conv: dict
        tuple -> int tax ids to integer representing their combo
    Returns
    -------
    None
    """

    for k, v in tax_id_2_prob.items():
        orm, created = EstimatedAbundance.objects.get_or_create(
            tax_id=k, task=task, name=names[k]
        )
        orm.abundance = v
        orm.save()


def save_uncertainty_prob(uncertainties, task, classification_converter):
    """
    Save the uncertainity probabilites into the database
    Parameters
    ----------
    uncertainties: np.ndarray
        The uncertainity probabilites, 2 x n_species, first col is lower CI, 2nd upper CI
    task: reads.models.JobMaster
        The task to be displayed
    classification_converter: dict
        tuple -> int, combination of tax ids to their representative classification number number

    Returns
    -------
    None
    """
    lower_bounds, upper_bounds = uncertainties
    int_to_taxid = {v: k for k, v in classification_converter.items()}
    total_uncerts = lower_bounds.size
    for i in range(total_uncerts):
        low_bound, up_bound = lower_bounds[i], upper_bounds[i]
        tax_ids = str(int_to_taxid[i])
        class_number = i
        UncertaintyProbability.objects.update_or_create(
            tax_ids=tax_ids,
            task=task,
            classification_number=class_number,
            defaults={"upper_ci_value": upper_bounds, "lower_ci_value": lower_bounds},
        )


def grammy_fit(
    classifications, init, draws, species_number, class_probs, tolerance=10e-10
):
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
    class_probs: np.ndarray
        the simulated probabilities that determine the probability a read of a certain
        species will receive a certain classification. Each row represents the true species of the read,
        and each column a classification where the column number corresponds to the classification number
        assigned according to the scheme above. Generated from centrifuge classificaiton of reads

    Returns
    -------
    return pi
    """
    err = 1
    start_prob = init
    z = np.zeros((len(classifications), species_number))
    while err > tolerance:
        for i, classification in enumerate(classifications):
            # class_prob = c_all[:, classification]
            kp = class_probs[:, classification]
            # todo check with michael what kr, kp, ki represent
            kr = kp * start_prob
            zp = kr / np.sum(kr)
            z[i,] = zp
        # @michael todo what do these names mean
        N = z.sum(axis=0)
        pil = N / draws
        err = np.sum(np.absolute(start_prob - pil))
        start_prob = pil
    return start_prob


def grammy_uncertainty(class_prob, prob_est, classifications):
    """

    Parameters
    ----------
    class_prob: np.ndarry
        probabilities of classification as each class for each species
    prob_est: np.ndarray
        estimated fragment proportions - output of grammy fitting procedure
    classifications: np.ndarray
        is the read data (so a vector with the classification of each read)

    Returns
    -------
    tuple
        output is a tuple first element is 95% lower bound, second element 95% upper bound
    """
    # todo what do all these represent?!
    lp = prob_est.size
    lphat = prob_est.size - 1
    phat = prob_est[:-1]
    K = len(classifications)
    IO = np.zeros((lphat, lphat))
    for i in range(lphat):
        for j in range(lphat):
            IOij = [0] * K
            for k in range(K):
                IOij1 = (
                    class_prob[i, classifications[k]]
                    - class_prob[lphat, classifications[k]]
                )
                IOij2 = (
                    class_prob[j, classifications[k]]
                    - class_prob[lphat, classifications[k]]
                )
                IOij3 = IOij1 * IOij2
                IOij4 = np.sum(phat * class_prob[:lphat, classifications[k]])
                IOij5 = (1 - np.sum(phat)) * class_prob[lphat, classifications[k]]
                IOij6 = (IOij4 + IOij5) ** 2
                IOij7 = IOij3 / IOij6
                IOij[k] = IOij7
            IO[i, j] = sum(IOij)
    IO = np.asfortranarray(IO)
    IOinv = lapack.dpotri(np.linalg.cholesky(IO).T)[0]
    est_se = np.sqrt(np.diag(IOinv))
    up_CI = phat + 1.96 * est_se
    low_CI = phat - 1.96 * est_se
    return (low_CI, up_CI)
