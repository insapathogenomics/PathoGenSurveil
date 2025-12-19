"""
This snakemake file comprises the set of rules which are necessary to run seqsero2 pipeline.

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

version = "1.0.0"
last_updated = "2025-09-01"

# 	settings	----------

working_dir = config["working_directory"]
fastq_flt_dir = working_dir + config["flt_fastq"]  # directory with trimmed files
sample_check_receipt = working_dir + "sample_check.txt"
seqsero_output_dir = working_dir + config["seqsero_output"]
seqsero_output = working_dir + config["seqsero_output"] + "SeqSero_result.tsv"


# 	rules	----------
rule all:
    input:
        seqsero_output,


rule seqsero2:
    input:
        fastq_flt_dir,
        sample_check_receipt,
    output:
        seqsero_output,
    conda:
        "/home/vpmixao/miniconda3/envs/seqsero2"
    message:
        "Running SeqSero2 pipeline..."
    resources:
        cpus_per_task=1,
        mem_mb=8000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "cat "
            + sample_check_receipt
            + " | while read id ; do SeqSero2_package.py -t 2 -m k -d "
            + seqsero_output_dir
            + "$id -i "
            + fastq_flt_dir
            + "$id\_R1\_001P.fastq.gz "
            + fastq_flt_dir
            + "$id\_R2\_001P.fastq.gz; done; grep -v 'Sample' "
            + seqsero_output_dir
            + "*/SeqSero_result.tsv > "
            + seqsero_output
        )
