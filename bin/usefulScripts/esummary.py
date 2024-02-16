"""
Calls Entrez.esummaty to build a DataFrame containing GSEs metadata
"""
from Bio import Entrez
import pandas as pd
import time
Entrez.email = "diegocar@lcg.unam.mx"
Entrez.api_key=  "34677fdcfd2f0659a7f9ee05ab6e44704f09"
from functions import format_time

def get_esumm_data(id, i, m,tries=10):
  print(f'Consulting id: {i}/{m}')
  try:
    with Entrez.esummary(db="gds", id=id, retmode='xml') as handle:
        esummary = Entrez.read(handle)[0]
        meta = {'id_num':id,
                'taxon': esummary['taxon'],
                'n_samples': int(esummary['n_samples']),
                'suppl_type':esummary['suppFile'],
                'PDAT':esummary['PDAT'],
                'GSM': ';'.join([f'3{"0"*(8-len(gsm:=i["Accession"].strip("GSM")))}{gsm}' for i in esummary['Samples']]),
                'GPL': ';'.join([f"1{'0'*(8-len(gpl))}{gpl}" for gpl in esummary['GPL'].split(";")])}
  except:
    if tries:
        return(get_esumm_data(id, tries-1))
    meta = {'id_num':id,
                'taxon': "Error",
                'n_samples': "Error",
                'suppl_type':"Error",
                'PDAT':"Error",
                'GSM' : "Error",
                'GPL':"Error"}
  finally: return(meta)

start = time.time()
with Entrez.esearch(db="gds", term="Expression profiling by array[GTYP] AND gse[ETYP]", retmax=100000) as handle:
   record = Entrez.read(handle)
id_list = record['IdList']

m = len(id_list)
esumm_data = [get_esumm_data(id,i,m) for i,id in enumerate(id_list)]

print("Building final table")
pd.DataFrame.from_dict(esumm_data).to_csv("../results/231018GEOmetadata/esummary_metadata_2.tsv", sep='\t')
print(f"Done in {format_time(time.time()-start)}")

#3985242