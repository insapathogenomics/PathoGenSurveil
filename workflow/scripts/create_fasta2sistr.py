#!/usr/bin/env python

"""
This script creates the tsv file necessary for sistr rule

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import pandas

version = "1.0.0"
last_updated = "2025-09-01"

def create_fasta2sistr(fasta, fasta2sistr):
    """ creates the fasta2sistr file """

    with open(fasta2sistr, "w+") as out:
        with open(fasta) as infile:
            lines = infile.readlines()
            for line in lines:
                l = line.split(".")[0]
                sample = l.split("/")[-1]
                print(str(sample) + "\t" + line.split("\n")[0], file = out)

if __name__ == "__main__":
    print("**********   running create_fasta2sistr.py  ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    create_fasta2sistr(snakemake.input[0], snakemake.output[0])