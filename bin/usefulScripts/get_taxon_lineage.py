"""
Calls Entrez.esummary to build a DataFrame containing the metadata for a list of GSMs accession numbers
"""
from Bio import Entrez
import pandas as pd
from time import time,sleep
Entrez.email = "diegocar@lcg.unam.mx"
Entrez.api_key=  "34677fdcfd2f0659a7f9ee05ab6e44704f09"
Entrez.sleep_between_tries = 0.5
from functions import format_time
import argparse
import logging
from random import choice
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)

def get_taxid_lineage(txid, i, m, tries = 10):
    logging.info("Consulting id %s/%s", i, m)
    try:
        with Entrez.efetch(db="taxonomy", id=txid) as handle:
            record = Entrez.read(handle)
        temp = {"TaxId": txid}
        temp.update({d["Rank"] : d["ScientificName"] for d in record[0]["LineageEx"]})
        return(temp)
    except:
        if tries: return(get_taxid_lineage(txid, i, m , tries -1))
        logging.warning("Error: %s", txid)
        return({})

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="Path to file with GSM list")
parser.add_argument("-o", "--out", help="File to print output")
args = parser.parse_args()

start = time()
with open(args.path) as file:
    txids = [txid.strip('\n') for txid in file.readlines()]

m = len(txids)
lineage = [get_taxid_lineage(txid, i+1,m) for i,txid in enumerate(txids)]

pd.DataFrame.from_dict(lineage).to_csv(args.out, sep='\t')
logging.info("Done in %s", format_time(time()-start))