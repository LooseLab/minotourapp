import pandas as pd
DATAFRAME_NAME = "metrics_df.p.gz"

def get_or_create_metrics_df(dir_to_write_to):
    """
    Get or create the metrics dataframe to be appended to.
    Parameters
    ----------
    dir_to_write_to: pathlib.PosixPath
        The directory that the artic results are being appended to.
    Returns
    -------

    """
    dataframe_destination = dir_to_write_to / DATAFRAME_NAME
    if not dataframe_destination.exists():
        df_old = pd.DataFrame()
    else:
        try:
            df_old = pd.read_pickle(str(dataframe_destination))
        except EOFError as e:
            df_old = pd.DataFrame()
    return df_old


def write_out_artic_metrics(df_new, df_old, dir_to_write_to):
    """
    Write out the artic metrics as a pickled dataframe
    Parameters
    ----------
    df_new: pandas.core.frame.DataFrame
        The dataframe of the metrics from the latest iteration, indexed by barcodes present
    df_old: pandas.core.frame.DataFrame
        The dataframe of the metrics from previous iterations, to be appended to.
    dir_to_write_to: pathlib.Posixpath
        The path to write the dataframe to

    Returns
    -------

    """
    if df_old.empty:
        df_old = df_new
    else:
        df_old = df_old.append(df_new)
    if not df_old.empty:
        df_old.to_pickle(str(dir_to_write_to / DATAFRAME_NAME))

