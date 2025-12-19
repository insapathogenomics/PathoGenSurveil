"""
This snakemake file comprises the set of rules which are necessary to run confindr pipeline.

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

version = "1.0.0"
last_updated = "2025-09-01"

# 	settings	----------

fastq_flt_dir = (
    config["working_directory"] + config["flt_fastq"]
)  # directory with trimmed files
confindr_report = (
    config["working_directory"] + config["confindr_output"] + "confindr_report.csv"
)
pass_fasta = config["working_directory"] + config["pass_fasta"]
checked_confindr_receipt = config["working_directory"] + "confindr_receipt.txt"


# 	rules	----------
rule all:
    input:
        checked_confindr_receipt,


rule confindr:
    input:
        fastq_flt_dir,
    output:
        confindr_report,
    conda:
        config["pathogensurveil_conda"]
    message:
        "Running confindr pipeline..."
    resources:
        cpus_per_task=4,
        mem_mb=16000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "confindr -i "
            + fastq_flt_dir
            + " -t {resources.cpus_per_task} -o "
            + config["working_directory"]
            + config["confindr_output"]
        )


rule update_pass_fasta:
    input:
        confindr_report,
        pass_fasta,
    output:
        checked_confindr_receipt,
    conda:
        config["pathogensurveil_conda"]
    message:
        "Updating assemblies passing QC..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/update_assemblies_pass.py"
