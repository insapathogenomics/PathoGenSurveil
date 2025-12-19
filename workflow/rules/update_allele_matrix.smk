"""
This snakemake file comprises the set of rules which are necessary to create or update the global allele matrix

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

version = "1.0.0"
last_updated = "2025-09-01"

# 	settings	----------

allele_matrix = config["alleles"]  # global allele matrix
run_allele_matrix = (
    config["working_directory"] + config["tag_general_output"] + "_alleles.tsv"
)  # final run allele matrix
new_allele_matrix = (
    config["working_directory"] + config["tag_general_output"] + "_global_alleles.tsv"
)


# 	rules	----------
rule all:
    input:
        new_allele_matrix,


rule create_new_allele_matrix:
    input:
        run_allele_matrix,
    output:
        new_allele_matrix,
    conda:
        config["pathogensurveil_conda"]
    message:
        "Creating the new allele matrix..."
    params:
        run_allele_matrix,
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/update_allele_matrix.py"
