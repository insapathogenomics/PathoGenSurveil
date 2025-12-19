#!/usr/bin/env python

"""
This script checks for duplicate samples and changes sample names.

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import sys
import glob
import pandas

version = "1.0.0"
last_updated = "2025-09-01"

# 	functions	----------

def check_samples(sample_info, fastq_dir, fastq1_extension, fastq2_extension, sample_check_receipt):
    new_name = {}
    fastq_dir = fastq_dir + "/"
    df = pandas.read_table(sample_info, dtype = str)
    duplicate_sample_ids = df[df["Sample_ID"].duplicated(keep=False)]["Sample_ID"].unique()
    if len(duplicate_sample_ids) > 0:
        print("Duplicate Sample_IDs found:", list(duplicate_sample_ids))
        sys.exit("Aborting PathoGenSurveil due to duplicate Sample_IDs.")
    else:
        print("No duplicate Sample_IDs found. Continuing PathoGenSurveil...")
    new_name = df.set_index("NGS_ID")["Sample_ID"].to_dict()
    
    samples_not_found = []
    out_file = open(sample_check_receipt, "w+")
    for filename in glob.glob(fastq_dir + "*" + fastq1_extension):
        fastq_name = filename.split(fastq_dir)[1]
        NGS_ID = fastq_name.split(fastq1_extension)[0]
        print(NGS_ID)
        if NGS_ID in new_name.keys():
            print(str(new_name[NGS_ID]), file = out_file)
            os.system("mv " + fastq_dir + NGS_ID + fastq1_extension + " " + fastq_dir + str(new_name[NGS_ID]) + "_R1_001.fastq.gz")
            os.system("mv " + fastq_dir + NGS_ID + fastq2_extension + " " + fastq_dir + str(new_name[NGS_ID]) + "_R2_001.fastq.gz")
        else:
            print("Sample " + str(NGS_ID) + " is not in the provided metadata! Please add it!")
            samples_not_found.append(NGS_ID)
    out_file.close()

    if len(samples_not_found) > 0:
        os.system("rm " + sample_check_receipt)

if __name__ == "__main__":
    print("**********   running check_samples.py  ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    check_samples(snakemake.input[0], snakemake.input[1], snakemake.params[0], snakemake.params[1], snakemake.output[0])
