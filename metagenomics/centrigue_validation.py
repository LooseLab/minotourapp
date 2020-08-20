import pandas as pd

# from alignment.tasks_alignment import align_reads
from metagenomics.models import (
    CentrifugeOutput,
    MappingTarget,
)


def separate_target_reads(df, task):
    """
    Separate out any reads from the dataframe that classified as targets
    Parameters
    ----------
    df: pandas.core.frame.DataFrame
         The centrifuge output dataframe
    task: reads.models.JobMaster
        The task ORM that is representing this metagenomics task
    Returns
    -------

    """
    flowcell = task.flowcell
    barcodes = df["barcode_name"].unique().tolist()
    # Add an all reads string to the list
    barcodes.append("All reads")
    targets = (
        MappingTarget.objects.filter(target_set=task.target_set)
        .values_list("species", flat=True)
        .distinct()
    )
    print("Flowcell id: {} - Targets are {}".format(task.flowcell.id, targets))
    # Get the target reads rows in a seperate dataframe
    name_df = df[df["name"].isin(targets)]
    previous_targets = CentrifugeOutput.objects.filter(
        species__in=targets, task=task, barcode_name__in=barcodes
    ).values("species", "tax_id", "num_matches", "barcode_name", "sum_unique")
    # Get the results into  a dataframe
    num_matches_per_target_df = pd.DataFrame.from_records(previous_targets.values())
    print(
        "Flowcell id: {} - The previous number of matches dataframe is {}".format(
            flowcell.id, num_matches_per_target_df
        )
    )