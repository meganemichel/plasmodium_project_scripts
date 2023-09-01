#!/usr/bin/env python
# coding: utf-8

import allel
import pandas as pd
import numpy as np
import zarr
import argparse
from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map



def read_zarr(zarr_path): 
    """Read .zarr files and output variant and callset tables."""
    callset = zarr.open_group(zarr_path, mode='r')
    variant_table = allel.VariantChunkedTable(callset['variants'], \
    names = ['CHROM', 'POS', 'REF', 'ALT', 'RegionType', 'FILTER_PASS', 'numalt', 'VQSLOD', 'CDS', 'is_snp'])
    return callset, variant_table

def filter_ind(meta_file, **kwargs): 
    """Read Pf6 sample provenance and sequencing metadata text file into DataFrame.
    Filter metadata DataFrame to exclude both strains failing QC and strains isolated from returning travelers."""
    filter_travelers = kwargs.get('filter_travelers', None)
    filter_populations = kwargs.get('filter_populations', None)
    samples = pd.read_csv(meta_file, sep = '\t')
    if filter_travelers: 
        sample_selection = (samples['QC pass'] == True) & (samples['Is returning traveller'] == False)
    else: 
        sample_selection = (samples['QC pass'] == True)

    if filter_populations: 
        sample_selection = (samples['Population'] )
    
    print ("{} individuals in the dataset pass QC and are included in subsequent analyses.".format(np.count_nonzero(sample_selection)))
    
    samples_subset = samples[sample_selection]
    samples_subset.reset_index(drop = True, inplace = True)
    return sample_selection, samples_subset

def variant_filter(variant_table, filter_expression): 
    """Filter variant table output by read_zarr fuction. 
    
    Function removes postions meeting any of the following criteria: 
    1. Positions annotated as failing QC (for details see https://www.malariagen.net/resource/26)
    2. Positions falling outside of coding sequences
    3. Multiallelic SNPs
    4. Positions with VQSLOD < 6
    """

    filter_pass = variant_table.eval(filter_expression)[:]
    
    ## Selection variants that satisfy both conditions
    variant_selection = filter_pass
    print ("{} variants pass filtering thresholds and are included in subsequent analyses.".format(np.count_nonzero(variant_selection)))
    variants_subset = variant_table.compress(variant_selection)
    return variant_selection, variants_subset

def vectorized_pseudohap(geno, seed): 
    """Takes 3-d numpy genotype array with reference and alternate allelic depths.
    Weights alleles according to number of observations and generates pseudohaploid genotype."""
    rng = np.random.default_rng(seed)
    # Reference allele encoded as 2, alternate allele encoded as 0
    ref = geno[:,:,0]
    alt = geno[:,:,1]
    null = (ref + alt == 0).astype(int)
    pseudohap = lambda r, a, n: rng.choice(r * ['2'] + a * ['0'] + n * ['9']).astype('int8')
    vfunc = np.vectorize(pseudohap)
    output = vfunc(ref, alt, null)
    return output

def parallelize_dataframe(array, func, n_cores, child_states):
    """Function taken from https://towardsdatascience.com/make-your-own-super-pandas-using-multiproc-1c04f41944a1.
    Paralellizes DataFrame operations across multicore machines.
    """
    #Split array into n parts
    array_split = np.array_split(array, n_cores)
    concat = np.concatenate(process_map(func, array_split, child_states, max_workers = n_cores, chunksize = 4))
    return concat

def exclude_uninformative(pseudo):
    """Takes 2-d pseudohaploidized array. 
    Removes positions that are not segregating in the population of interest. 
    Also removes singletons that are uninformative for genetic analysis."""
    reference = np.count_nonzero(pseudo == 2, axis=1)
    alternate = np.count_nonzero(pseudo == 0, axis=1)
    uncalled = np.count_nonzero(pseudo == 9, axis=1)
    ac = np.column_stack((reference, alternate, uncalled))
    bool_array = ac[:, :2].min(axis=1) > 1

    print ("{} segregating variants are non-singletons and are included in subsequent analyses.".format(np.count_nonzero(bool_array)))

    return pseudo[bool_array, :], bool_array
    
def np2eigenstrat(geno_array, samples_subset, variants_subset, out_path, base, snp_bool_array):
    """Outputs genetic data in eigenstrat format for downstream analysis.
    
    Takes the following inputs: 
    1. geno_array- pseudoahploid numpy array output by pseudohap function
    2. samples_subset- samples_subset dataframe output by the filter_ind function
    3. variants_subset- .zarr compressed table containing filtered variants from variant_filter function
    4. out_path- path for output eigenstrat data
    5. base- str for eigenstrat file names
    """
    ##Convert genotype array to eigenstrat .geno file
    pseudo = geno_array
    np.savetxt(out_path + base + '.geno',            pseudo , fmt='%s', delimiter = '')  
    
    ##Convert population metadata array to eigenstrat .ind file
    df = pd.DataFrame(samples_subset)[['Sample', 'Population']]
    df.insert(loc = 1, column = 'Sex', value = 'U')
    df.to_csv(out_path + base + '.ind', sep = '\t', index = False, header = False)
    
    ##Extract SNP data from variant table and write to eigenstrat .snp file
    chrom = variants_subset['CHROM'][:]
    pos = variants_subset['POS'][:]
    ref = variants_subset['REF'][:]
    alt = [item[0] for item in list(variants_subset['ALT'])]
    temp_snp_data = pd.DataFrame([chrom, pos, ref, alt], index = ['CHROM', 'POS', 'REF', 'ALT']).transpose()
    s = pd.Series(list(snp_bool_array), name='bools')
    snp_data = temp_snp_data[s.values]
    snp_data.insert(loc = 1, column = 'GENETIC_POS', value = '0')
    snp_data.to_csv(out_path + base + '.snp', sep = '\t', index = True, header = False)
    
def np2eigenstrat_nosubset(geno_array, samples_subset, variants_subset, out_path, base):
    """Outputs genetic data in eigenstrat format for downstream analysis.
    
    Takes the following inputs: 
    1. geno_array- pseudoahploid numpy array output by pseudohap function
    2. samples_subset- samples_subset dataframe output by the filter_ind function
    3. variants_subset- .zarr compressed table containing filtered variants from variant_filter function
    4. out_path- path for output eigenstrat data
    5. base- str for eigenstrat file names
    """
    ##Convert genotype array to eigenstrat .geno file
    pseudo = geno_array
    np.savetxt(out_path + base + '.geno', pseudo, fmt='%s', delimiter = '')  
    
    ##Convert population metadata array to eigenstrat .ind file
    df = pd.DataFrame(samples_subset)[['Sample', 'Population']]
    df.insert(loc = 1, column = 'Sex', value = 'U')
    df.to_csv(out_path + base + '.ind', sep = '\t', index = False, header = False)
    
    ##Extract SNP data from variant table and write to eigenstrat .snp file
    chrom = variants_subset['CHROM'][:]
    pos = variants_subset['POS'][:]
    ref = variants_subset['REF'][:]
    alt = [item[0] for item in list(variants_subset['ALT'])]
    snp_data = pd.DataFrame([chrom, pos, ref, alt], index = ['CHROM', 'POS', 'REF', 'ALT']).transpose()
    snp_data.insert(loc = 1, column = 'GENETIC_POS', value = '0')
    snp_data.to_csv(out_path + base + '.snp', sep = '\t', index = True, header = False)

def rename_chromosomes(out_path, base): 
    """Reformats chromosome names to comply with smartPCA format specifications."""
    new_name = str(out_path + base) + '.reformat.snp'
    df = pd.read_csv(out_path + base + '.snp', sep = '\t', names = ['SNP', 'Chromosome', 'Genetic_Distance', 'Physical_Distance', 'Reference', 'Alternate'])
    chromosomes = {'Pf3D7_01_v3':'1', 'Pf3D7_02_v3': '2', 'Pf3D7_03_v3': '3', 'Pf3D7_04_v3': '4', 'Pf3D7_05_v3': '5', \
     'Pf3D7_06_v3': '6', 'Pf3D7_07_v3':'7', 'Pf3D7_08_v3': '8', 'Pf3D7_09_v3': '9', 'Pf3D7_10_v3': '10', \
     'Pf3D7_11_v3': '11', 'Pf3D7_12_v3': '12', 'Pf3D7_13_v3': '13', 'Pf3D7_14_v3': '14', 'PvP01_01_v1': '1', \
    'PvP01_02_v1': '2', 'PvP01_03_v1': '3', 'PvP01_04_v1': '4', 'PvP01_05_v1': '5', 'PvP01_06_v1': '6', \
    'PvP01_07_v1': '7', 'PvP01_08_v1': '8', 'PvP01_09_v1': '9', 'PvP01_10_v1': '10', 'PvP01_11_v1': '11', \
    'PvP01_12_v1': '12', 'PvP01_13_v1': '13', 'PvP01_14_v1': '14'}
    df.insert(loc = 1, column = 'Reformat_Chrom', value = [chromosomes[item] for item in list(df.Chromosome)])
    df.drop( labels = ['Chromosome'], axis = 1, inplace = True)
    df.to_csv(new_name, sep = '\t', index = False, header = False)

def write_bed(outpath, base): 
    """Rewrites .snp file according to bedfile specifications."""
    new_name = str(outpath + base) + '.bed'
    df = pd.read_csv(outpath + base + '.reformat.snp', sep = '\t', names = ['SNP', 'Chromosome', 'Genetic_Distance', 'Physical_Distance', 'Reference', 'Alternate'])
    chromosome = list(df['Chromosome'])
    start = [int(pos) - 1 for pos in list(df['Physical_Distance'])]
    stop = [int(pos) for pos in list(df['Physical_Distance'])]
    test = ['test'] * df['Chromosome'].shape[0]
    qual = ['1000'] * df['Chromosome'].shape[0]
    plus = ['+'] * df['Chromosome'].shape[0]
    df = pd.DataFrame([chromosome, start, stop, test, qual, plus], index = ['Chromosome', 'Start', 'Stop', 'Test', 'Qual', 'Strand'],\
                  columns = None).transpose()
    df.to_csv(new_name, sep='\t', index = False, header = None)


def main(): 
    parser = argparse.ArgumentParser(description = "Parameters for ZarrToEigenstrat.")
    parser.add_argument('zarr_dir', help = "Directory containing Zarr files.")
    parser.add_argument('eigenstrat_outdir', help = "Directory for eigenstrat format output files.")
    parser.add_argument('base', help = "")
    parser.add_argument('metadata_file', help = "Path to metadata file.")
    parser.add_argument('filter_expression', help = "Filter expression for pruning variant set.")
    parser.add_argument('-f', '--filtertravelers', help = "Remove strains isolated from returning travelers.", action=argparse.BooleanOptionalAction)
    parser.add_argument('-s', '--snpsubset', help = "Remove uninformative positions, including both SNPs not segregating in samples of interest and singleton positions.", action=argparse.BooleanOptionalAction)
    parser.add_argument('-n', '--numpy_seed', help = "Numpy seed for pseudohaplotyping.", required = False, default = "")
    parser.add_argument('-c', '--cores', help = "Cores available for parallelization steps.", required = False, default = 1)
    try: 
        args = parser.parse_args()
    except: 
        parser.print_help()
        sys.exit(0)

    if args.numpy_seed: 
        print ('Setting random number generator with seed {}.'.format(args.numpy_seed))
        rng = np.random.default_rng(int(args.numpy_seed))

    else: 
        print ('Default random number generator used.')
        rng = np.random.default_rng(None)

    ss = rng.bit_generator._seed_seq
    child_states = ss.spawn(int(args.cores))

    print ('Reading Zarr files.')
    ##Read .zarr data
    callset, variant_table = read_zarr(args.zarr_dir)

    print ('Filtering variants.')
    ##Filter callset variants to include low quality variant calls
    variant_selection, variants_subset = variant_filter(variant_table, args.filter_expression)
    
    print ('Removing individuals failing QC.')
    ##Fitler callset samples to exclude samples failing QC and isolates from recent travelers
    if args.filtertravelers: 
        sample_selection, samples_subset = filter_ind(args.metadata_file, filter_travelers = True)
    else: 
        sample_selection, samples_subset = filter_ind(args.metadata_file)
    
    print ('Subsetting genotypes based on previous filtering steps.')
    ##Subset genotype data based on previous filtering steps
    genotypes = allel.GenotypeChunkedArray(callset['calldata']['AD'])
    geno = genotypes.subset( variant_selection, sample_selection )
    geno = geno[:, :, :2].astype('int16')

    print ('Generating pseudohaploid genotypes.')
    pseudohap_df = parallelize_dataframe(geno, vectorized_pseudohap, int(args.cores), child_states)

    if args.snpsubset: 
        print ('Check and remove uninformative positions.')
        pseudo_subset, bool_array = exclude_uninformative(pseudohap_df)
        np2eigenstrat(pseudo_subset, samples_subset, variants_subset, args.eigenstrat_outdir, args.base, bool_array)
    else: 
        np2eigenstrat_nosubset(pseudohap_df, samples_subset, variants_subset, args.eigenstrat_outdir, args.base)

    print ('Renaming chromsomes.')
    #Rename chromosomes to conform to eigenstrat format
    rename_chromosomes(args.eigenstrat_outdir, args.base)
    write_bed(args.eigenstrat_outdir, args.base)

if __name__ == "__main__": 
    main()

