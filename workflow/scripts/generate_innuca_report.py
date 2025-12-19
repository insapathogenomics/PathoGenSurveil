#!/usr/bin/env python

"""
This script summarizes and reports the INNUca run

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import glob
import pandas

version = "1.0.1"
last_updated = "2025-09-16"

def join_df(old_df,new_df):
	""" join two dataframes
	input: pandas dataframe
	output: pandas dataframe
	"""
	
	old_df.set_index(old_df.columns[0], inplace = True)
	new_df.set_index(new_df.columns[0], inplace = True)
	final_df = pandas.concat([old_df, new_df])
	final_df = final_df.reset_index()

	return final_df

def preparing_report(scripts_folder, samples_report, final_report, fasta2allele):
    """ summarize and report INNUca run """

    final_samples_report = pandas.DataFrame()
    final_combined_report = pandas.DataFrame()

    for report in samples_report:
        sample_tag = report.split("_samples_report.tsv")[0]
        samples_report = pandas.read_table(report, dtype = str)
        combined_report = pandas.read_table("".join([sample_tag, "_combine_samples_report.tsv"]), dtype = str)
        samples_report = samples_report.set_index(samples_report.columns[0], drop = True)
        combined_report = combined_report.set_index(combined_report.columns[0], drop = False)
        if final_samples_report.empty:
            final_samples_report = samples_report
        else:
            final_samples_report = join_df(final_samples_report, samples_report)
        if final_combined_report.empty:
            final_combined_report = combined_report
        else:
            final_combined_report = join_df(final_combined_report, combined_report)
        #os.system("rm " + report + " " + sample_tag + "_combine_samples_report.tsv")
    final_report_df = pandas.concat([final_combined_report, final_samples_report], axis=1)

    columns_report = []
    with open(scripts_folder + "/columns_innuca_report.txt") as infile:
        lines = infile.readlines()
        for line in lines:
            l = line.split("\n")[0]
            columns_report.append(l)

    final_report_df = final_report_df[columns_report]
    final_report_df.to_csv(final_report, index = False, header = True, sep = "\t")

    pass_samples = final_report_df[final_report_df["samples_passQC"] != "FAIL"]
    final_assemblies = pass_samples["final_assembly"].values.tolist()

    with open(fasta2allele, "w+") as outfile:
        for fasta in final_assemblies:
            print(fasta, file = outfile)

if __name__ == "__main__":
    print("**********   running generate_innuca_report.py  ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    preparing_report(snakemake.input[0], snakemake.input[1:], snakemake.output[0], snakemake.output[1])
