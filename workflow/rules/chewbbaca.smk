"""
This snakemake file comprises the set of rules which are necessary to run INSA's allele calling pipeline for foodborne bacteria.
This pipeline uses chewBBACA and requires its docker installation. https://github.com/B-UMMI/chewBBACA

By Veronica Mixao
@INSA
vmixao@gmail.com
"""

version = "1.0.0"
last_updated = "2025-09-01"


# 	settings	----------

files4allele_calling = (
    config["working_directory"] + config["pass_fasta"]
)  # file with assemblies for chewBBACA
clean_allele_matrix = (
    config["working_directory"] + config["tag_general_output"] + "_alleles.tsv"
)  # final allele matrix


# 	functions	----------


def params_allecall(dictionary):
    """This function obtains all the parameters necessary to run chewBBACA from the config file"""

    command_line_chewbbaca = ""
    for param in config["chewbbaca"]:  # getting chewBBACA parameters
        if (
            config["chewbbaca"][param] != "None"
            and config["chewbbaca"][param] != False
            and config["chewbbaca"][param] != "off"
        ):
            if config["chewbbaca"][param] == True:
                command_line_chewbbaca += " --" + param
            else:
                command_line_chewbbaca += (
                    " --" + param + " " + str(config["chewbbaca"][param])
                )

    return command_line_chewbbaca


# 	rules chewBBACA	----------


rule all:
    input:
        clean_allele_matrix,


rule chewbbaca:
    input:
        files4allele_calling,
    output:
        config["working_directory"] + config["allele_output"] + "results_alleles.tsv",
    conda:
        config["pathogensurveil_conda"]
    message:
        "Running chewBBACA pipeline..."
    resources:
        cpus_per_task=config["chewbbaca"]["cpu"],
        mem_mb=int(config["chewbbaca"]["cpu"]) * 4000,
        slurm_extra=config["slurm_extra"],
    shell:
        (
            "singularity exec -B "
            + config["mount"]
            + " "
            + config["chewbbaca_container"]
            + " chewie AlleleCall -i {input[0]} -g "
            + config["chewbbaca_schema"]
            + " -o "
            + config["working_directory"]
            + config["allele_output"]
            + params_allecall(config)
            + "; mv "
            + config["working_directory"]
            + config["allele_output"]
            + "*/results_alleles.tsv "
            + config["working_directory"]
            + config["allele_output"]
        )


rule clean_chewie_matrix:
    input:
        config["working_directory"] + config["allele_output"] + "results_alleles.tsv",
    output:
        clean_allele_matrix,
    # conda:
    #    config["pathogensurveil_conda"]
    message:
        "Cleaning the allele matrix..."
    resources:
        cpus_per_task=1,
        mem_mb=4000,
        slurm_extra=config["slurm_extra"],
    script:
        "../scripts/clean_chewie_matrix.py"
