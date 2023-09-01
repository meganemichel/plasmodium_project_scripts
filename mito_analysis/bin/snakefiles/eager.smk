################################################################################
# Project: Plasmodium Nextstrain
# Title: Eager module
#
# Desciption: Run downloaded files through nf-core/eager to generate SNP alignments.
#
# Megan Michel, 08/05/2021
################################################################################


print (scripts_dir)
rule run_nfcore:
    """
    Pre-process and map fastq data to a reference genome with nf-core/eager.
    """
    params: 
        config = config['CONFIG'],
        profiles = config['PROFILES'],
        email = config['EMAIL'], 
        vcfs = config['ADDITIONAL_VCFs'],
        regions_exclude = config['REGIONS_EXCLUDE'],
        reference_gff = config['REFERENCE_GFF'], 
        vcf_dir = config['SSLIB_vcfs'], 
        genoSL_table = config['GENOSL_TABLE']
    input: 
        ["input_files/genome_to_read/{sample}_reads.fastq.gz".format(sample=sample) for sample in SAMPLES]
    output: 
        snpTable = "eager/multivcfanalyzer/snpTable.tsv.gz",
        snpAlignment = "eager/multivcfanalyzer/snpAlignment.fasta.gz"
    run:
        in_list = expand("input_files/genome_to_read/{sample}_reads.fastq.gz", sample = SAMPLES)
        shell("mkdir -p ./vcfs")
        shell('for item in $(find {params.vcfs} -name "*vcf.gz"); do echo $item; ln -sf $item ./vcfs; done')
        if config['SSLIB_vcfs']: 
            shell('for item in $(find {params.vcf_dir} -name "*vcf.gz"); do ln -sf $item ./vcfs; done')
        shell("nextflow run nf-core/eager -r 2.4.6 -c {params.config} -profile {params.profiles} \
             --outdir ./eager -work-dir ./work --input './input_files/genome_to_read/*.fastq.gz' \
             --run_multivcfanalyzer --min_base_coverage 3 --additional_vcf_files './vcfs/*vcf.gz' \
              --reference_gff_annotations {params.reference_gff} \
              --reference_gff_exclude {params.regions_exclude} --single_end \
            -with-tower --email {params.email}  --skip_qualimap")
        shell("mv .nextflow* eager")
        shell("mv work eager")
        
rule gather_vcfs:
    """
    Pre-process and map fastq data to a reference genome with nf-core/eager.
    """
    params: 
        genoSL_table = config['GENOSL_TABLE']
    input:
        snp_table = rules.run_nfcore.output.snpTable, 
        snp_alignment = rules.run_nfcore.output.snpAlignment
    output: 
        snpTable = "multivcf_out/snpTable.tsv",
        snpAlignment = "multivcf_out/snpAlignment.fasta"
    run: 
        shell("mkdir -p multivcf_out")
        if config['SSLIB_vcfs']:
            shell("zcat eager/multivcfanalyzer/snpTable.tsv.gz > eager/multivcfanalyzer/snpTable.tsv")
            shell("{scripts_dir}/GenoSL/genoSL.R ./eager/multivcfanalyzer/snpTable.tsv {params.genoSL_table} --output ./multivcf_out/snpTable")
            shell("mv ./multivcf_out/snpTable_genotyped.fasta ./multivcf_out/snpAlignment.fasta")
            shell("mv ./multivcf_out/snpTable_genotyped.tsv ./multivcf_out/snpTable.tsv")
        else: 
            shell("cp eager/multivcfanalyzer/snpTable.tsv.gz multivcf_out")
            shell("gzip -d multivcf_out/snpTable.tsv.gz")
            shell("cp eager/multivcfanalyzer/snpAlignment.fasta.gz multivcf_out")
            shell("gzip -d multivcf_out/snpAlignment.fasta")
