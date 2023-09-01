import pandas as pd
import csv
import numpy as np

def filter_tsv(file, bad_samples, out): 
    # Read tsv file into a dataframe
    df = pd.read_csv(file, sep = '\t')
    df.columns = [col.replace('_reads.unifiedgenotyper.vcf','') for col in list(df.columns)]
    df.columns = [col.replace('.unifiedgenotyper.vcf','') for col in list(df.columns)]

    print ('columns', df.columns)
    # Read tsv file with samples to exclude into list
    with open(bad_samples, newline='') as f: 
        reader = csv.reader(f)
        data = list (reader)
    data = [item.strip() for sublist in data for item in sublist]
    data = [item.replace('.unifiedgenotyper.vcf','') for item in data]
    data = [item.replace('_reads', '') for item in data]
    print ('bad sample names', data)
    
    # Subset dataframe to exclude samples in list
    df.drop(labels = data, axis = 1, inplace = True)
        
    #Write dataframe to new .tsv file
    df.to_csv(out, sep='\t', header = True, index = False)

def snps_only(multi_vcf_snp_table): 
    """Read SNP table into a DataFrame and filter positions for which only one allele is called."""
    snp_table = pd.read_csv(multi_vcf_snp_table, sep = '\t', index_col = 0, header = 0)
    temp = snp_table.replace("N", np.nan)
    # Get list of positions where there are at least two different allele calls
    uniqueValues = temp.nunique()
    uniqueValues = uniqueValues[uniqueValues > 1]
    good_pos = list(uniqueValues.index)
    # Subset SNP table to include only positions in good_pos list
    snp_table = snp_table[good_pos]
    return snp_table

def write_fasta(snp_table, fasta_out): 
    """Takes dataframe from snps_only function and outputs fasta file"""
    df = snp_table

    with open(fasta_out, 'w') as f: 
        for index, row in df.iterrows(): 
            print (index)
            f.write('>' + index)
            f.write('\n')
            seq = ''.join(list(row))
            f.write(seq)
            f.write('\n')





    