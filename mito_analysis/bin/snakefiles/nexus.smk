################################################################################
# Project: Plasmodium Nextstrain
# Title: Nexus format module
#
# Desciption: 
#
# Megan Michel, 08/05/2021
################################################################################
import nexus_format
import snp_table_filter
import os
import pandas as pd

rule format_fasta:
    """
    Reformat multifasta file to exclude reference sequence.
    """
    input: 
        snpTable =  rules.gather_vcfs.output.snpTable,
        snpAlignmet = rules.gather_vcfs.output.snpAlignment
    output:
        snpAlignment_noReference = temp("alignments/fullAlignment_noReference.fasta")
    shell:
        """
        mkdir -p ./temp_fastas
        faSplit byname {input.snpAlignmet} temp_fastas/
        cat $(find temp_fastas -name "*fa" | grep -v Reference) > {output.snpAlignment_noReference}
        rm -r temp_fastas
        """

rule low_cov_calc: 
    """
    Identify individuals with genome coverage lower than specified threshold.
    """
    params: 
        missingness_threshold = config['MISSINGNESS_THRESHOLD'],
        outgroup = config['OUTGROUP']
    input: 
        snpAlignment = rules.format_fasta.output.snpAlignment_noReference
    output: 
        low_cov_samples = temp("low_cov_samples.txt"),
        nuc_freq = temp("nuc_freq.txt")
    run: 
        shell("faCount {input} > {output.nuc_freq}")
        df = pd.read_csv(output[1], sep = '\t')
        df['Missing_Percent'] = [(ns / length)* 100 for ns, length in zip(list(df.N), list(df.len))]
        print (df)
        print(float(params.missingness_threshold))
        sub = df[df['Missing_Percent'] > float(params.missingness_threshold)] 
        print (len(sub))
        with open(output.low_cov_samples, 'w') as f: 
            for index, row in sub.iterrows():
                if index == df.index[-1]:
                    print('This row is last')
                else:
                    print (row["#seq"])
                    print (str(row["#seq"].replace('.unifiedgenotyper.vcf', '').replace('_reads', '')))
                    f.write(str(row["#seq"].replace('.unifiedgenotyper.vcf', '').replace('_reads', '')))
                    f.write("\n")
        f.close() 

rule remove_bad_samples: 
    """
    Reformat SNP table to exclude low coverage individuals. 
    """
    input: 
        low_cov_samples = rules.low_cov_calc.output.low_cov_samples, 
        snpTable = rules.gather_vcfs.output.snpTable
    output: 
        filtered_tsv = "alignments/snpTable/snpTable_exclude_lowcov_indivs.tsv"
    run:
        snp_table_filter.filter_tsv(input[1], input[0], output[0])

rule alignment_filter: 
    params: 
        partial_deletion_filter = config['PARTIAL_DELETION_FILTER'],
        exclude = config['TO_EXCLUDE']
    input: 
        filtered_tsv = rules.remove_bad_samples.output.filtered_tsv
    output: 
        nexdir = directory("alignments/partial_deletion_filter")
    shell:
        """
        mkdir {output.nexdir}
        if test -f "{params.exclude}"; then
            echo "{params.exclude} exists."
            Rscript {scripts_dir}/MDF/MDF.R -t {params.exclude} {input.filtered_tsv} {params.partial_deletion_filter} filtered  > alignments/partial_deletion_filter/filtered.log
        else
            echo "No {params.exclude}."
            Rscript {scripts_dir}/MDF/MDF.R {input.filtered_tsv} {params.partial_deletion_filter} filtered  > alignments/partial_deletion_filter/filtered.log
        fi
        mv filtered* alignments/partial_deletion_filter
        """

rule snp_only: 
    params: 
        partial_deletion_filter = config['PARTIAL_DELETION_FILTER']
    input: 
        nexdir = rules.alignment_filter.output.nexdir
    output: 
        finaldir = directory("alignments/final"), 
        temp_tsv = temp('tmp.tsv')
    run: 
        fasta = shell('cp $(find alignments/partial_deletion_filter -name "*tsv") tmp.tsv')
        df = snp_table_filter.snps_only(output[1])
        snp_num = str(df.shape[1])
        partial_deletion_filter = str(params.partial_deletion_filter)
        shell("mkdir alignments/final")
        df.to_csv("alignments/final/filtered_" + partial_deletion_filter + "_pd_" + snp_num + "_positions.snps_only.tsv", sep='\t', header = True, index = False)
        snp_table_filter.write_fasta(df, "alignments/final/filtered_" + partial_deletion_filter + "_pd_" + snp_num + "_positions.snps_only.fasta")


rule make_nexus: 
    input: 
        finaldir = rules.snp_only.output.finaldir
    output: 
        filtered_nex = temp("filtered_snp_alignment.nex")
    shell:
        """
        file=$(find {input.finaldir} -name "*fasta")
        seqmagick convert --output-format nexus --alphabet dna $file {output.filtered_nex}
        """

rule add_trait_data: 
    params: 
        fixed_meta = config['FIXED_META']
    input: 
       filtered_nex = rules.make_nexus.output.filtered_nex, 
    output: 
        final_nex = "nexus_final.nex", 
    run:
        nexus_format.make_nexus(input[0], params.fixed_meta)
        os.rename(input[0], output[0])



