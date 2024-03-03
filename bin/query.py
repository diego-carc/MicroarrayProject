'''
Gets GEO metadata
'''
from Bio import Entrez
import re
import pandas as pd
from multiprocessing import Pool
import time
import argparse
import logging
import os
from utils_query import *
from sys import stdout


### MAIN ###
# Log config
total_time = time.time()
os.system("export GEOPARSE_USE_HTTP_FOR_FTP=yes")
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO, stream=stdout)
sep = \
"""
##############################################################################
{}
##############################################################################
"""

# Argumentos
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxID", help="path of file with taxID list")
parser.add_argument("-o", "--outFile", help="name of output file")
parser.add_argument("-p", "--path", help="path to output dir", default='./')
parser.add_argument("-e", "--email", help="email to use with Bio.Entrez")
parser.add_argument("-a", "--apiKey", help="If needed, the api Key to use in Bio.Entrez", required=False)
parser.add_argument("-s", "--sleepTime", help="Sleep between tries used in Bio.Entrez run", default=1, type=int)
parser.add_argument("-c", "--cores", help="Cores to be used", default=1, type=int)
parser.add_argument("-m", "--max", help="Max number of Entrez retrieve", default=100000, type=int)
parser.add_argument("-r", "--rank", help="Specificity rank level", default="sameTaxId", choices=["specie", "sameTaxId"])
args = parser.parse_args()

# Asignación de argumentos
logging.info(sep.format("Start running"))
arguments = \
  f"""Run with the following arguments:
    --taxID: {args.taxID}
    --outFile: {(output := args.outFile)}
    --path: {(path := args.path)}
    --email: {args.email}
    --apiKey: {(args.apiKey)}
    --sleep_between_tries: {args.sleepTime}
    --cores: {args.cores}
    --max: {args.max}
    --rank: {args.rank}
"""
with open(args.taxID) as file:
  txid_list = [txid.rstrip('\n') for txid in file.readlines()]
Entrez.email = args.email
Entrez.api_key = args.apiKey
Entrez.sleep_between_tries = args.sleepTime
logging.info("Argument reading succesful")
logging.info(arguments)

# Obtener los GEO accession
logging.info(sep.format("Start GEO consulting"))
start = time.time()
GSE_IDs = [(id,re.sub(r'^20+', 'GSE', gse),path) for id in txid_list for gse in get_GSE_IDs(id, args.max)]

logging.info("Total consulting time: %s", format_time(time.time()-start))

# Obtener las tablas de metadatos por txID
logging.info(sep.format("Start metadata parsing"))
start = time.time()
with Pool(args.cores) as p:
  metadata_sync = p.map(get_GEO_metadata_from_txid, GSE_IDs)

logging.info("Total metadata parsing time: %s", format_time(time.time()-start))

raw_metadata_table = pd.concat(metadata_sync).reset_index().drop('index', axis=1)
raw_metadata_table.GPL_manufacturer = raw_metadata_table.GPL_manufacturer.apply(rename_manufacturer)
raw_metadata_table.to_csv(os.path.join(path, f"raw_{output}"), sep= '\t', header= True, index= False)

# Filtrar contaminaciones
logging.info(sep.format("Filtering"))
all_metadata_table = apply_GEO_filters(raw_metadata_table, args.rank, path, output)
logging.info("Filtering done")

# Stats
logging.info("Building output...")
all_metadata_table['TaxId'].value_counts().plot(kind='bar', legend=True).get_figure().savefig(os.path.join(path, "GSMs_txID.png"))

for name in ['GPL_technology', 'GPL_distribution', 'GSM_type', 'GSM_channel_count']:
  get_plot(all_metadata_table, name, path)

rel = pd.read_csv(os.path.join(path, f"taxonomyRelation_{output}"), sep='\t')

# Data type
all_metadata_table.TaxId = all_metadata_table.TaxId.astype("int64")
raw_metadata_table.TaxId = raw_metadata_table.TaxId.astype("int64")
rel.TaxId = rel.TaxId.astype("int64")

plot_species_dist(all_metadata_table, rel).get_figure().savefig(os.path.join(path, "GSMs_species_distribution.png"))

try:
  comp_filt(all_metadata_table, "Filtered", raw_metadata_table, "Raw", rel).get_figure().savefig(os.path.join(path, "raw_filtered_comp.png"))
except: logging.warning("Unable to plot raw_filtered_comp.png")

# Construir archivo output
all_metadata_table.to_csv(os.path.join(path, output), sep= '\t', header= True, index= False)
logging.info("Total time: %s", format_time(time.time()-total_time))
logging.info("Thanks for using")

#apikey: 34677fdcfd2f0659a7f9ee05ab6e44704f09