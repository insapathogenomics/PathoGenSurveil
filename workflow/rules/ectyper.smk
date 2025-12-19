"""
This snakemake file comprises the set of rules which are necessary to run ECTyper pipeline.

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

version = "1.0.0"
last_updated = "2025-09-30"

# 	settings	----------

fasta_files = config["working_directory"] + config["pass_fasta"]
ectyper_report = config["working_directory"] + config["ectyper_output"] + "output.tsv"
ectyper_final_report = (
    config["working_directory"] + config["tag_general_output"] + "_ectyper.tsv"
)
assembly_directory = config["working_directory"] + "assemblies/"


# 	rules	----------
rule all:
    input:
        ectyper_report,
        ectyper_final_report,


rule create_assembly_directory:
    input:
        fasta_files,
    output:
        directory(assembly_directory),
    conda:
        config["pathogensurveil_conda"]
    message:
        "Creating assembly directory..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "mkdir "
            + assembly_directory
            + "; cat "
            + fasta_files
            + " | while read fasta ; do cp $fasta "
            + assembly_directory
            + " ; done"
        )


rule ectyper:
    input:
        directory(assembly_directory),
    output:
        ectyper_report,
    conda:
        config["pathogensurveil_conda"]
    message:
        "Running ectyper pipeline..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "ectyper -i "
            + assembly_directory
            + " -o "
            + config["working_directory"]
            + config["ectyper_output"]
            + " --pathotype --verify"
        )


rule clean_ectyper:
    input:
        ectyper_report,
    output:
        ectyper_final_report,
    # conda:
    #    config["pathogensurveil_conda"]
    message:
        "Cleaning the Ectyper matrix..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/clean_ectyper.py"
