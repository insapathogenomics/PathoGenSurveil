"""
This snakemake file comprises the set of rules which are necessary to run sistr pipeline.

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

version = "1.0.0"
last_updated = "2025-09-01"

# 	settings	----------

fasta_files = config["working_directory"] + config["pass_fasta"]
sistr_directory = config["working_directory"] + config["sistr_output"]
sistr_report = config["working_directory"] + config["tag_general_output"] + "_sistr.tsv"
fasta_sistr = config["working_directory"] + "fasta2sistr.tsv"


# 	rules	----------
rule all:
    input:
        sistr_report,


rule create_fasta2sistr:
    input:
        fasta_files,
    output:
        fasta_sistr,
    conda:
        config["pathogensurveil_conda"]
    message:
        "Creating input for sistr pipeline..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/create_fasta2sistr.py"


rule sistr:
    input:
        fasta_sistr,
    output:
        sistr_report,
    conda:
        config["pathogensurveil_conda"]
    message:
        "Running sistr pipeline..."
    resources:
        cpus_per_task=1,
        mem_mb=8000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "mkdir "
            + sistr_directory
            + "; cat "
            + fasta_sistr
            + " | while read id fasta ; do sistr --qc -vv --alleles-output "
            + sistr_directory
            + "$id\_allele-output.csv --cgmlst-profiles "
            + sistr_directory
            + "$id\_cgmlst-profiles.csv -f tab -o "
            + sistr_directory
            + "$id\_sistr-output.tab $fasta; done; grep -v 'cgmlst_ST' "
            + sistr_directory
            + "*sistr-output.tab > "
            + sistr_report
        )
