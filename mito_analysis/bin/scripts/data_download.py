from Bio import SeqIO
from Bio import Entrez
import pandas as pd

def ncbi_search_nucleotide(search):
    """Search ncbi nucleotide database for records matching user-supplied search term."""
    handle = Entrez.esearch(db = 'nucleotide', idtype = 'acc', retmax= 100000, term= search)
    records = Entrez.read(handle)
    handle.close()
    return records

def ncbi_download(IdList, gb_file): 
    """Write output of ncbi_search_nucleotide to genbank file."""
    gb_hd = Entrez.efetch(db = 'nucleotide', id = IdList, rettype = 'gb')
    data = gb_hd.read()
    out_handle = open(gb_file, "w")
    out_handle.write(data)
    out_handle.close()

def get_metadata(gb_file): 
    """Parse genbank file and extract relevant metadata."""
    input_handle  = open(gb_file, "r") 
    data = []
    col_names = ['Id', 'Name', 'Database References', 'Description', 'Molecule Type', 'Topology', 'Data File Division', \
            'Date', 'Accessions', 'Sequence Version', 'Keywords', 'Source', 'Organism', 'Taxonomy', 'References', \
            'Organelle', 'Molecule Type', 'Isolate', 'Isolation Source', 'Host', 'db_xref', 'Country', 'lat_lon', \
            'Collection Date', 'Collected By']
    for index, record in enumerate(SeqIO.parse(input_handle, 'genbank')): 
        genome = [record.id, record.name, '; '.join(record.dbxrefs), record.description, \
                   record.annotations['molecule_type'], record.annotations['topology'], \
                   record.annotations['data_file_division'], record.annotations['date'], \
                   '; '.join(record.annotations['accessions']), record.annotations['sequence_version'], 
                  '; '.join(record.annotations['keywords']), record.annotations['source'], \
                   record.annotations['organism'], '; '.join(record.annotations['taxonomy'])]
        references = []
        for item in record.annotations['references']: 
                ref = item.authors + ' ' + item.title + '. ' + item.journal + '.'
                references.append(ref)
        genome.append('; '.join(references))
        for feature in record.features: 
            if feature.type == 'source': 
                features = [(feature.qualifiers.get('organelle', ''), feature.qualifiers.get('mol_type', ''), \
                      feature.qualifiers.get('isolate', ''), feature.qualifiers.get('isolation_source', ''), \
                      feature.qualifiers.get('host', ''), feature.qualifiers.get('db_xref', ''), \
                       feature.qualifiers.get('country', ''), feature.qualifiers.get('lat_lon', ''), \
                    feature.qualifiers.get('collection_date', ''), feature.qualifiers.get('collected_by', ''))]
                for sample in features: 
                    feature_anno = ['; '.join(item) for item in sample ]
                genome = genome + feature_anno
        data.append(genome)
    df = pd.DataFrame(data, columns = col_names)
    return df

def fasta_download(IdList, fasta_file): 
    """Download multifasta file containing sequences for list of accession numbers."""
    fasta = Entrez.efetch(db = 'nucleotide', id = IdList, rettype = 'fasta')
    data_fasta = fasta.read()
    out_fasta = open(fasta_file, "w")
    out_fasta.write(data_fasta)
    out_fasta.close()