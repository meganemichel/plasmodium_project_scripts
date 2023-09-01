################################################################################
# Project: Plasmodium Nextstrain
# Title: Data download module
#
# Description: Downloads fasta data from the NCBI Nucleotide database matching
# given search term. Parses metadata from genbank files and returns pseudogenized
# reads from downloaded data.
#
# Megan Michel, 08/05/2021
################################################################################

import data_download
import pandas as pd
from Bio import Entrez

search = config['SEARCH_TERM']
Entrez.email =  config['ENTREZ_EMAIL']

#Identify records matching search term for download
records = data_download.ncbi_search_nucleotide(search)['IdList']
SAMPLES = records + [config['OUTGROUP']]
print (len(SAMPLES))

rule ncbi_download: 
    """
    Download fasta and metadata files from NCBI.
    """
    output:
        gb = temp("temp.gb"),
        fa = "concat.fa"
    run:
        data_download.ncbi_download(SAMPLES, output.gb)
        data_download.fasta_download(SAMPLES, output.fa)
        
rule compile_metadata: 
    """
    Reformat metadata file.    
    """
    input: 
        gb =  rules.ncbi_download.output.gb
    output: 
        meta = "input_files/ncbi_download_meta.csv"
    run: 
        df = data_download.get_metadata(input.gb)
        with open(output.meta, 'w') as f: 
            df.to_csv(f, index = False)
        f.close()

rule fasta_split: 
    """
    Split multifasta by sample.
    """
    input: 
        fasta = rules.ncbi_download.output.fa
    output: 
        expand("input_files/sample_fastas/{sample}.fa", sample=SAMPLES)
    shell: 
        """
        faSplit byname {input.fasta} input_files/sample_fastas/
        mv concat.fa input_files
        """

rule genome_to_read:
    """
    Generate simulated reads from fasta files.
    """
    input: 
        fa = "input_files/sample_fastas/{sample}.fa"
    output: 
        fq = "input_files/genome_to_read/{sample}_reads.fastq.gz"
    shell: 
        """
        Genome2Reads.jar {input.fa} 50 1 | gzip > {output.fq}
        """



    

    