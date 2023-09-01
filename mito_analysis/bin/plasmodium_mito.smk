################################################################################
# Project: Plasmodium Mitochondrial Analysis
#
#	NOTE: THIS IS COMMAND LINE FORMAT: 
#	snakemake --configfile ./config.json --cores 4 -s plasmodium_mito.smk 
#
#
# Megan Michel, 08/05/2021
################################################################################

import os
import data_download

# -----------------------------------------------------------------------------#
#                                 Setup                                        #
# -----------------------------------------------------------------------------#
workdir: config['RESULTS_DIR']

if not config:
    print("ERROR: Please specify --configfile")
    quit(1)

pipeline_dir = os.path.dirname(workflow.basedir)
rules_dir = os.path.join(pipeline_dir,"bin", "snakefiles")
scripts_dir = os.path.join(pipeline_dir, "bin", "scripts")

#Sub snakefiles
include: rules_dir + "/data_download.smk",
include: rules_dir + "/eager.smk", 
include: rules_dir + "/nexus.smk"

# -----------------------------------------------------------------------------#
# Main Target                                                                  #
# -----------------------------------------------------------------------------#

rule all:
    """
    The default pipeline targets.
    """
    input: 
        # expand("input_files/genome_to_read/{sample}_reads.fastq.gz", sample=SAMPLES), 
        #"input_files/ncbi_download_meta.csv", 
        # expand("eager/genotyping/{sample}_reads.unifiedgenotyper.vcf.gz", sample=SAMPLES),
        # "eager/multivcfanalyzer/fullAlignment.fasta.gz", 
        # #"filtered_snp_alignment.nex"
        "nexus_final.nex" 

