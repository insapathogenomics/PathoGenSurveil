"""
This snakemake file comprises the set of rules which are necessary to run INSA's genome assembly pipeline for foodborne bacteria.
This pipeline uses INNUca and requires its docker installation. https://github.com/B-UMMI/INNUca

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

import glob

version = "1.0.0"
last_updated = "2025-09-01"


# 	settings	----------

working_dir = config["working_directory"]
fastq_dir = working_dir + config["fastq"]  # directory with fastq files
innuca_report = (
    working_dir + config["tag_general_output"] + "_innuca_report.tsv"
)  # final report
files4allele_calling = (
    working_dir + config["pass_fasta"]
)  # file with assemblies for chewBBACA
fastq_flt_dir = working_dir + config["flt_fastq"]  # directory with trimmed files
sample_check_receipt = working_dir + "sample_check.txt"


#   variables   ----------

SAMPLES = []
with open(sample_check_receipt) as insamples:
    lines = insamples.readlines()
    for line in lines:
        l = line.split("\n")[0]
        SAMPLES.append(l)


# 	functions	----------


def params_innuca(dictionary):
    """determine INNUca's command line"""

    command_line_innuca = ""
    for param in config["innuca"]:  # getting INNUca parameters
        if (
            config["innuca"][param] != "None"
            and config["innuca"][param] != False
            and config["innuca"][param] != "off"
        ):
            if config["innuca"][param] == True:
                command_line_innuca += " --" + param
            else:
                command_line_innuca += (
                    " --" + param + " " + str(config["innuca"][param])
                )
    return command_line_innuca


# 	rules INNUCA	----------


rule all:
    input:
        innuca_report,
        files4allele_calling,
        directory(fastq_flt_dir),


rule innuca:
    input:
        fastq_dir + "{sample}_R1_001.fastq.gz",
        fastq_dir + "{sample}_R2_001.fastq.gz",
    output:
        config["working_directory"]
        + config["assembly_output"]
        + "{sample}_samples_report.tsv",
        config["working_directory"]
        + config["assembly_output"]
        + "{sample}_combine_samples_report.tsv",
    conda:
        config["pathogensurveil_conda"]
    message:
        "Running INNUca pipeline on sample {wildcards.sample}..."
    resources:
        cpus_per_task=config["innuca"]["threads"],
        mem_mb=config["innuca_mem_mb"],
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "singularity exec -B "
            + config["mount"]
            + " "
            + config["innuca_container"]
            + ' python /NGStools/INNUca/INNUca.py -f {input[0]} {input[1]} -s "'
            + config["species"]
            + '" -o '
            + config["working_directory"]
            + config["assembly_output"]
            + "{wildcards.sample}/"
            + params_innuca(config)
            + "; mv "
            + config["working_directory"]
            + config["assembly_output"]
            + "{wildcards.sample}/samples_report.*.tab "
            + config["working_directory"]
            + config["assembly_output"]
            + "{wildcards.sample}_samples_report.tsv; mv "
            + config["working_directory"]
            + config["assembly_output"]
            + "{wildcards.sample}/combine_samples_reports.*.tab "
            + config["working_directory"]
            + config["assembly_output"]
            + "{wildcards.sample}_combine_samples_report.tsv"
        )


rule innuca_report:
    input:
        config["scripts"],
        expand(
            config["working_directory"]
            + config["assembly_output"]
            + "{sample}_samples_report.tsv",
            sample=SAMPLES,
        ),
    output:
        innuca_report,
        files4allele_calling,
    # conda:
    #    config["pathogensurveil_conda"]
    message:
        "Getting INNUca reports..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/generate_innuca_report.py"


rule create_flt_fastq:
    input:
        innuca_report,
    output:
        directory(fastq_flt_dir),
    # conda:
    #    config["pathogensurveil_conda"]
    message:
        "Generating directory with filtered fastq..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "mkdir "
            + fastq_flt_dir
            + "; mv "
            + config["working_directory"]
            + config["assembly_output"]
            + "*/*/trimmomatic/*_001P.fastq.gz "
            + fastq_flt_dir
        )
