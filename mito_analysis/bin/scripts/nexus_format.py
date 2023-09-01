#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import dendropy
import numpy as np
import pycountry
import pycountry_convert as pc

def make_nexus(infile, metadata_file): 
    ## VARIABLE DEFINITION
    file = open(infile)

    ## Extract sample names in order from Nexus file
    nex = file.read()
    d = dendropy.DataSet()
    d.read_from_string(nex, "nexus")
    taxa = d.taxon_namespaces[0]
    
    ## Read metadata file into dataframe
    df = pd.read_csv(metadata_file, sep =',', index_col = 'Id')
    df['Country'].replace(np.nan, 'unknown', inplace = True)

    ## Get geographic region for each country in dataframe
    geo = {}
    for item in set([item.split(':')[0] for item in df.Country]): 
        try: 
            a2 = pycountry.countries.search_fuzzy(item)[0].alpha_2
            geo[item] = pc.country_alpha2_to_continent_code(a2)
        except: 
            if item == 'South_America': 
                geo[item] = 'SA'
            elif item == 'Central_America': 
                geo[item] = 'NA'
            elif item == 'Asia': 
                geo[item] = "AS"
            elif item == 'S.C.America': 
                geo[item] = 'unknown'
            else: 
                geo[item] = item
    
    ## Add region data to dataframe and reformat
    region = [geo[country] for country in [item.split(':')[0] for item in df.Country]]
    new_df = pd.DataFrame([list(df.index), list(df['Country']), region], index = ['Name', 'Country', 'Region']).transpose()
    
    ## Reorder regions to match sample order in nexus file
    geo = []
    for sample in taxa: 
        taxon = str(sample).strip("'").replace(' ', '_')
        geo.append(new_df[new_df.Name == taxon].Region.item())
    
    ## Generate text to write to output file
    ntraits = len(set(geo))
    traits = list(set(geo))
    trait_labels = '\tTraitLabels'
    for label in traits: 
        label = label.replace(' ', '_')
        trait_labels = trait_labels + ' ' + str(label)
    trait_labels = trait_labels + ';'
    to_write = ['BEGIN TRAITS;', '\tDimensions NTRAITS=' + str(ntraits) +';', \
        '\tFormat labels=yes missing=? separator=Spaces;', trait_labels, '\tMatrix']
    for taxon, zone in zip(taxa, geo):
        string = str(taxon).strip("'").replace(' ', '_')
        for label in traits: 
            if zone == label: 
                string = string + ' 1'
            else: 
                string = string + ' 0'
        to_write.append(string)
    to_write.append(';')
    to_write.append('END;')
    
    ## Write trait data to nexus file
    with open(infile, 'a') as f: 
        for line in to_write: 
            f.write(line)
            f.write('\n')
    f.close()