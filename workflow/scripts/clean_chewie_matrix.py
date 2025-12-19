#!/usr/bin/env python

"""
This script cleans the allele matrix provided by chewBBACA

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import glob
import pandas

version = "1.0.0"
last_updated = "2025-09-01"

def clean_matrix(original, final):
    """ cleans the allele matrix provided by chewBBACA """

    missing_data = {"PLOT3": "0", "PLOT5": "0", "LOTSC": "0", "NIPH": "0", "NIPHEM": "0", "ALM": "0", "ASM": "0", "PAMA": "0", "LNF": "0"}

    mx = pandas.read_table(original, dtype=str)
    mx[mx.columns[0]][0] = mx[mx.columns[0]][0].split(".contigs.length_GCcontent_kmerCov.mappingCov.excluded_contigs.fasta")[0]
    mx = mx.applymap(lambda x: missing_data.get(x, x))
    mx = mx.replace("INF-","", regex=True)
    mx = mx.replace("\*","", regex=True)

    mx.to_csv(final, index = False, header = True, sep = "\t")

if __name__ == "__main__":
    print("**********   running clean_chewie_matrix.py  ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    clean_matrix(snakemake.input[0], snakemake.output[0])