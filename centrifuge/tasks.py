# # TODO maybe rewrite the centrifuger task into this file
# from celery import task
# from centrifuge.models import CentrifugeOutput, CentrifugeOutputParsed
# from jobs.models import JobMaster
# import pandas as pd
# import numpy as np
#
#
# def create_centrifuge_models(row, classified_per_barcode):
#     """
#     Append a CentOutput object to a list for each row in the centrifuge output
#     :param row: the row from the data frame
#     :param classified_per_barcode: The number of reads classified for each barcode as a dictionary
#     :return: The list of newly created objects
#
#     """
#     if row["proportion_of_classified"] == "Unclassified" or row["proportion_of_classified"] == np.nan:
#         row["proportion_of_classified"] = round(row["num_matches"] / classified_per_barcode[row["barcode_name"]], 3)
#     # logger.info(row["proportion_of_classified"])
#     # logger.info(row["barcode_name"])
#     # logger.info(row["proportion_of_classified"] == np.NaN)
#     # logger.info(row["proportion_of_classified"] == np.nan)
#     return CentrifugeOutput(name=row["name"],
#                                   tax_id=row["tax_id"],
#                                   task=row["task"],
#                                   num_matches=row["num_matches"],
#                                   sum_unique=row["sum_unique"],
#                                   barcode_name=row["barcode_name"],
#                                   proportion_of_classified=row["proportion_of_classified"],
#                                   superkingdom=row["superkingdom"],
#                                   phylum=row["phylum"],
#                                   classy=row["class"],
#                                   order=row["order"],
#                                   family=row["family"],
#                                   genus=row["genus"],
#                                   species=row["species"],
#                                   latest=row["latest"])
#
#
# @task()
# def output_parser(task_id):
#     """
#     Parse the centrifuge output data and save into a site readable format
#     :param task_id: The primary key for this task on this flowcell
#     :return:
#     """
#
#     task = JobMaster.objects.get(pk=task_id)
#
#     flowcell = task.flowcell
#
#     metagenomics_task = JobMaster.objects.get(flowcell=flowcell, job_type__name="Metagenomics")
#
#     new_data = CentrifugeOutput.objects.filter(process=False, task=metagenomics_task)
#
#     parsed_data = CentrifugeOutputParsed.objects.filter(task=metagenomics_task)
#
#     parsed_data_df = pd.DataFrame(list(parsed_data.values()))
#
#     new_data_df = pd.DataFrame(list(new_data.values()))
#
#     new_data_gb = new_data_df.groupby(["name", "barcode_name"])
#
#     new_data_df.set_index(["name", "barcode_name"], inplace=True)
#
#     new_data_df["num_matches"] = new_data_gb["num_matches"].sum()
#
#     new_data_df["sum_unique"] = new_data_gb["sum_unique"].sum()
#
#     new_data_df.reset_index(inplace=True)
#
#     new_data_df.drop_duplicates(subset=["name", "barcode_name"], inplace=True)
#
#     if not parsed_data_df.empty:
#         new_latest = parsed_data_df["latest"].unique().tolist()[0] + 1
#
#         merged_data_df = pd.merge(parsed_data_df, new_data_df, on=["barcode_name", "name", "superkingdom",
#                                                                    "phylum", "classy", "order", "family",
#                                                                    "genus", "species", "tax_id"],
#                                   how="outer")
#
#         values = {"num_matches_x": 0, "num_matches_y": 0, "sum_unique_x": 0, "sum_unique_y": 0}
#
#         merged_data_df.fillna(value=values, inplace=True)
#
#         merged_data_df["task"] = metagenomics_task
#
#         merged_data_df["num_matches"] = merged_data_df["num_matches_x"] + merged_data_df["num_matches_y"]
#
#         merged_data_df["sum_unique"] = merged_data_df["sum_unique_x"] + merged_data_df["sum_unique_y"]
#
#         gb_bc = merged_data_df.groupby("barcode_name")
#
#         classed_per_barcode = gb_bc["num_matches"].sum().to_dict()
#
#         barcodes = merged_data_df["barcode_name"].unique().tolist()
#
#         to_save_df = merged_data_df
#     else:
#         new_latest = 1
#
#         gb_bc = new_data_df.groupby("barcode_name")
#
#         classed_per_barcode = gb_bc["num_matches"].sum().to_dict()
#
#         barcodes = new_data_df["barcode_name"].unique().tolist()
#
#         to_save_df = new_data_df
#
#     def divd(row, cl_bar):
#         """
#         Calculate_proportion_of_classified
#         :param row: The df row
#         :param cl_bar: The number of reads classified in a barcode, in dict form keyed to the barcode_name
#         :return:
#         """
#         return round((row["num_matches"] / cl_bar[row["barcode_name"]]) * 100, 4)
#
#     to_save_df["proportion_of_classified"] = to_save_df.apply(divd, args=(classed_per_barcode,), axis=1)
#
#     to_save_df["latest"] = new_latest
#
#     to_save_df["task"] = metagenomics_task
#
#     to_save_df.rename(columns={"classy": "class"}, inplace=True)
#
#     to_bulk_save_series = to_save_df.apply(create_centrifuge_models, args=(classed_per_barcode,), axis=1)
#
#     CentrifugeOutputParsed.objects.bulk_create(to_bulk_save_series.values.tolist(), batch_size=20000)
#
#
#     metagenomics_task.latest_batch = new_latest
#
#     metagenomics_task.save()
#
#     CentrifugeOutputParsed.objects.filter(task=metagenomics_task, latest__lt=new_latest).delete()
#
#
#
