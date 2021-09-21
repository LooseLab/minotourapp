"""utils.py for the useful logic functions for readfish views"""
import numpy as np
import pandas as pd


def get_start_to_start_coords(jm):
    """
    Get the starting coordinates for this primer scheme
    Parameters
    ----------
    jm: reads.models.JobMaster
        JobMaster instance for artic
    Returns
    -------
    pd.core.frame.DataFrame
        Dataframe of ranges for both forward and reverse strand
    """
    bed = pd.read_csv(
        jm.primer_scheme.bed_file,
        sep="\t",
        header=None,
        names=["chromosome", "start", "stop", "primer", "group"],
    )
    bed[["amplicon", "strand"]] = bed["primer"].str.split("_", expand=True)[[1, 2]]
    bed["strand"] = bed["strand"].map({"LEFT": "+", "RIGHT": "-"})
    bed["amplicon"] = bed["amplicon"].astype(int)
    bed["mod_start"] = bed.groupby("strand")["start"].transform(lambda g: g.shift(-1))
    bed["mod_stop"] = bed.groupby("strand")["stop"].transform(lambda g: g.shift())
    bed["range"] = np.where(bed["strand"] == "+", bed["mod_start"], bed["mod_stop"])
    bed["range_start"] = np.where(bed["strand"] == "+", bed["start"], bed["range"])
    bed["range_end"] = np.where(bed["strand"] == "+", bed["range"], bed["stop"])
    bed["range_start"] = bed["range_start"].fillna(0)
    bed["range_end"] = bed["range_end"].fillna(
        jm.reference.length
    )  # set to max length of reference
    return bed
