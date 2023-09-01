#!/usr/bin/env python
# coding: utf-8

# **Author:** Megan Michel <br>
# **Date:** 02/12/2021 <br>
# **Description:** Script used to filter output of MALT in gzipped Blast-TAB format to identify taxa that may cause false-positives during pathogen screening. <br>
# **Usage:** db_false_pos_taxa_id_021221.py <blast_tab.txt> <outfile.csv> < percent ID filter> <br>
# **Notes:** Unlike in previous iterations, if a read maps to two references from the same taxid, I keep only the accession with the highest percent identity. I also have a new method for removing alignments to Plasmodium spp. Together these measures should reduce the size of the Plasmodium screening database. 

# In[51]:


import gzip
import pandas as pd
import csv
from ete3 import NCBITaxa
ncbi = NCBITaxa()
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
import sys
import time
import csv
from Bio import Entrez
Entrez.email = "megan_michel@g.harvard.edu"
##Database updated last on 02/12/2021, base conda env
#ncbi.update_taxonomy_database()


# In[78]:


def read_df(file):
    """Take BLAST-Tab format file and read into pandas dataframe"""
    d = '/'
    name = file.split('/')[:-1]
    outdir = d.join(name)
    bad_taxids = get_ipython().getoutput('zcat $file | cut -f2 | grep -v "|" | sort |uniq ')
    print (bad_taxids)
    taxid_dict = {}
    bad_tax_str = ', '.join(bad_taxids)
    handle =  Entrez.esummary(db="nuccore", id=bad_tax_str, retmode="xml")
    records = Entrez.parse(handle)
    for acc, record in zip(bad_taxids, records):
    # each record is a Python dictionary or list.
        taxid = int(record['TaxId'])
        taxid_dict[acc] = taxid
    with gzip.open(file, "rt") as f:
        df = pd.DataFrame([column for column in csv.reader(f,delimiter='\t')], columns= ['Query_ID',     'Reference_ID','Percent_Identity', 'Alignmnet_Length', 'Mismatch_Number', 'Gap_Open',     'Alignment_Start_Query', 'Alignment_End_Query', 'Alignment_Start_Subject', 'Alignment_End_Subject',     'Evalue', 'Bitscore'])
    df['Percent_Identity'] = df['Percent_Identity'].astype(float)
    return df, taxid_dict

def percentID_filter(df, percent_ID): 
    """Take df from read_df, filter to contain entries above a certain percent identity."""
    df = df[(df['Percent_Identity']) > int(percent_ID)]
    return df

def reformat(df, taxid_dict): 
    """Take df from read_df and reformat to explicitly define Taxid and Taxonomy for references."""
    taxids = []
    taxa = []
    for index, row in df.iterrows(): 
        try: 
            taxid = int(row['Reference_ID'].split('|')[2])
            taxids.append(int(taxid))
            taxa.append(row['Reference_ID'].split('|')[0])            
        except IndexError: 
            taxid = taxid_dict[str(row['Reference_ID'])]
            taxids.append(taxid)
            taxa.append(str(row['Reference_ID']))
    df.insert(1, "Reference_Taxid", taxids, True) 
    df.insert(2, "Reference_Taxonomy", taxa, True) 
    df = df.drop("Reference_ID",  axis=1)
    return df

def filter_plasmodium(df): 
    """Takes deduplicated df from remove_duplicates. Removes taxids that are children to Plasmodium node."""
    descendants = ncbi.get_descendant_taxa(5820, intermediate_nodes = True)
    descendants.append(5820)
    for item in descendants: 
        df = df[~df['Reference_Taxid'].isin(descendants)]
    return df

def remove_duplicates(df): 
    """Take reformatted df and filter to keep only highest percent identity match for each
    Query/Taxid combination."""
    df.sort_values(by=['Query_ID', 'Reference_Taxid', 'Percent_Identity'], inplace=True, ascending=False)
    df.drop_duplicates(subset=['Query_ID', 'Reference_Taxid'], keep='first', inplace=True)
    return df
    
def summary_table(df): 
    """Takes deduplicated df from remove_duplicates. Generates file summarizing matches by taxon."""
    summary_list = []
    read_num = len(df.Query_ID.unique())
    alignment_num = df.shape[0]
    for taxid in df.Reference_Taxid.unique(): 
        sub = df[df['Reference_Taxid'] == taxid]
        taxon = ncbi.get_taxid_translator([taxid])
        percent = (int(sub.shape[0]) / int(alignment_num)) * 100
        max_pident = sub.Percent_Identity.max()
        mean_pident = sub.Percent_Identity.mean()
        AC_numbers = set(list(sub.Reference_Taxonomy))
        summary_list.append([taxid, taxon.get(int(taxid)), max_pident, mean_pident, ';'.join(AC_numbers)])
    df_fin = pd.DataFrame(summary_list, columns=['Taxid', 'Taxonomy', 'Maximum_Percent_ID', 'Mean_Percent_ID', 'Accession_Numbers'])
    return df_fin

def write_summary_table(df, outfile): 
    df.to_csv(outfile, index=False)
    


# In[95]:


# NOTE: CAN BE USED FOR TESTING/TROUBLESHOOTING
#df, taxid_dict = read_df('/projects1/users/michel/plasmodium/mirror_data/tab/test.fastq.gz')
# filter_df = percentID_filter(df, 90)
# reformat_df = reformat(filter_df, taxid_dict)
# plasmo_no_df = filter_plasmodium(reformat_df)
# dedup_df = remove_duplicates(reformat_df)
# final_df = summary_table(dedup_df)


# In[ ]:


def main(): 
    if len(sys.argv) < 3: 
        raise Exception("Usage: db_false_pos_taxa_id_021221.py <blast_tab.txt.gz> <outfile.csv> <percent ID filter>")
    infile, output_file, percent_id = sys.argv[1:]
    
    #Read BLAST-tab format file into dataframe. 
    df, taxid_dict = read_df(infile)
    print ('Loaded {} alignments from BLAST-tab alignment file output by MALT.'.format(df.shape[0]))
    print ('Obtained {} taxonic IDs for references lacking this infomration in the full NT database.'.format(len(taxid_dict)))

    #Filter to include only alignments with greater than 90% identity
    aln_df = percentID_filter (df, percent_id)
    print ('Retained {} alignments with greater than 90% identity.'.format(aln_df.shape[0]))

    #Reformat dataframe, creating separate fields for taxid and accession ID. 
    reformat_df = reformat(aln_df, taxid_dict)
    print ('Reformatted dataframe containing {} alignments'.format(reformat_df.shape[0]))
    
    #Filter dataframe to remove any alignments to taxa mapping to the Plasmodium node or its descendants.
    filter_df = filter_plasmodium(reformat_df)
    print ('Retained {} alignments from taxa other than Plasmodium spp.'.format(filter_df.shape[0]))
    
    #For each read/taxid combo, keep only best-matching reference
    dedup_df = remove_duplicates(filter_df)
    print ('Retained {} alignments with unique read and taxid combinations.'.format(dedup_df.shape[0]))

    #Summarized stats for each potential false positive taxon
    data = summary_table(dedup_df)
    print ('Generated report on alignments to {} potential false-positive taxa.'.format(data.shape[0]))
    print ('{} unique false-positive references will be included in final database.'.format(len([item for sublist in [item.split(';') for item in list(data.Accession_Numbers)] for item in sublist])))

    #Write summarized table to .csv file
    write_summary_table(data, output_file)


# In[ ]:


if __name__ == "__main__":
    main()

