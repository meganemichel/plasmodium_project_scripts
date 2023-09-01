#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import argparse

def read_eigenstrat(in_dir, in_prefix): 
    """Read .ind, .snp, and .geno files in eigenstrat format. Outputs pandas dataframes for each file."""
    ind = pd.read_csv(in_dir + in_prefix + '.ind', delim_whitespace=True, names = ['Ind', 'Sex', 'Pop'])
    snp = pd.read_csv(in_dir + in_prefix + '.snp', delim_whitespace=True,\
                      names = ['Ignore', 'Chromosome', 'Ignore_2', 'Position', 'Ref', 'Alt'])
    snp.index = [str(chr) + '_' + str(pos) for chr, pos in zip(snp.Chromosome, snp.Position)]
    geno_temp = pd.read_csv(in_dir + in_prefix + '.geno', names = ['temp'])
    geno = geno_temp['temp'].apply(lambda x: pd.Series(list(x)))
    geno.columns = list(ind.Ind)
    geno.index = [str(chr) + '_' + str(pos) for chr, pos in zip(snp.Chromosome, snp.Position)]
    return geno, snp, ind

def fix_missing(geno): 
    """Replace missing data with 'N' in genotype dataframe."""
    geno_out = geno.replace('9', 'N')
    return geno_out

def replace_ref_alt(geno, snp, ind):
    """Replace numerical SNP encodings with alleles based on SNP data."""
    rows = []
    for index, row in geno.iterrows(): 
        ref = snp.loc[index] ['Ref']
        alt = snp.loc[index] ['Alt']
        row = row.replace('0', alt)
        row = row.replace('2', ref)
        rows.append(list(row))
    geno_final = pd.DataFrame(rows, columns = ind.Ind)
    return geno_final

def eigen_to_fasta(geno, out_dir, out_prefix):   
    """Write output to multifasta format."""    
    with open(out_dir+ '/' + out_prefix + '.fasta', 'w') as f: 
        for ind in geno.columns: 
            f.write(">" + ind)
            f.write("\n")
            f.write(''.join(list(geno[ind])))
            f.write("\n")
    f.close()

def main(): 
    parser = argparse.ArgumentParser(description = "Parameters for eigen_to_fasta.")
    parser.add_argument('in_dir', help = "Directory containing eigenstrat files.")
    parser.add_argument('in_prefix', help = "Prefix for eigenstrat format files with endings .geno, .ind, and .snp.")
    parser.add_argument('out_dir', help = "Directory containing output files in multifasta format.")
    parser.add_argument('out_prefix', help = "Prefix for multifasta output file with ending .fasta.")
    args = parser.parse_args()
    
    geno, snp, ind = read_eigenstrat(args.in_dir , args.in_prefix)
    print (args.in_dir + args.in_prefix + '.geno')
    geno_out = fix_missing(geno)
    final = replace_ref_alt(geno_out, snp, ind)
    eigen_to_fasta(final, args.out_dir, args.out_prefix)
    
if __name__ == "__main__": 
    main()

