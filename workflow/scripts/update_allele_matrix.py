#!/usr/bin/env python

"""
This script creates or updates the allele matrix used for genomic surveillance

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import glob
import pandas

version = "1.0.0"
last_updated = "2025-09-01"

def update_allele_matrix(original, run, final):
    """ updates the allele matrix used for genomic surveillance """

    run_mx = pandas.read_table(run, dtype = str)
    if original == "":
        run_mx.to_csv(final, index = False, header = True, sep = "\t")
    else:
        mx_original = pandas.read_table(original, dtype = str)
        mx_final = pandas.concat([mx_original, run_mx], ignore_index=True)
        mx_final.to_csv(final, index = False, header = True, sep = "\t")

if __name__ == "__main__":
    print("**********   running update_allele_matrix.py  ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    update_allele_matrix(snakemake.params[0], snakemake.input[0], snakemake.output[0])