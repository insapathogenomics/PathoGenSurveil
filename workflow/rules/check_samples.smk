"""
This snakemake file comprises the set of rules which are necessary to check for duplicate samples and change sample names.

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


# 	settings	----------

working_dir = config["working_directory"]
sample_info = config["metadata"]
fastq_dir = working_dir + config["fastq"]  # directory with fastq files
sample_check_receipt = working_dir + "sample_check.txt"
fastq1_extension = config["fastq1_extension"]
fastq2_extension = config["fastq2_extension"]


# 	rules	----------


rule all:
    input:
        sample_check_receipt,


rule check_samples:
    input:
        sample_info,
        fastq_dir,
    output:
        sample_check_receipt,
    # conda:
    #    config["pathogensurveil_conda"]
    message:
        "Checking samples..."
    params:
        fastq1_extension,
        fastq2_extension,
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/check_samples.py"
