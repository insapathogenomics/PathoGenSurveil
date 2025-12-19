#!/usr/bin/env python

"""
This script cleans the fasta2chewbbaca file

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import os
import pandas

version = "1.0.0"
last_updated = "2025-09-01"

def clean_assemblies_file(report, fasta_file, receipt):
    """ cleans the fasta2chewbbaca file """

    mx = pandas.read_csv(report, dtype=str)
    sample_to_status = dict(zip(mx["Sample"], mx["ContamStatus"]))

    with open("temp_fasta2chewbbaca.txt", "w+") as out:
        with open(receipt, "w+") as receipt_file:
            with open(fasta_file) as fasta:
                lines = fasta.readlines()
                for line in lines:
                    l = line.split(".")[0]
                    sample = l.split("/")[-1]
                    if sample_to_status[sample] == "False":
                        print(line.split("\n")[0], file = out)
                    print(sample, sample_to_status[sample], file = receipt_file)
    
    os.system("mv temp_fasta2chewbbaca.txt " + fasta_file)

if __name__ == "__main__":
    print("**********   running update_assemblies_pass.py  ****")
    print("version:", version)
    print("last_updated", last_updated, "\n")
    clean_assemblies_file(snakemake.input[0], snakemake.input[1], snakemake.output[0])