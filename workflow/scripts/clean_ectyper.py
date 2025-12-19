#!/usr/bin/env python

"""
This script cleans the Ectyper matrix

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import glob
import pandas

version = "1.0.0"
last_updated = "2025-09-30"

def clean_matrix(original, final):
    """ cleans the ectyper matrix """
    
    mx = pandas.read_table(original, dtype=str)
    mx[mx.columns[0]] = mx[mx.columns[0]].str.split('.').str[0]
    mx = mx[~mx[mx.columns[0]].isna() & (mx[mx.columns[0]] != "")]

    mx.to_csv(final, index = False, header = True, sep = "\t")

if __name__ == "__main__":
    print("**********   running clean_ectyper.py   ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    clean_matrix(snakemake.input[0], snakemake.output[0])