"""
Calls Entrez.esummary to build a DataFrame containing the metadata for a list of GSMs accession numbers
"""
from Bio import Entrez
import pandas as pd
from time import time,sleep
Entrez.email = "diegocar@lcg.unam.mx"
Entrez.api_key=  "34677fdcfd2f0659a7f9ee05ab6e44704f09"
from functions import format_time
import argparse
import logging
from random import choice
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)

def get_sum(id, i,m, tries=10):
    sleep(choice([0.1,0.5,1]))
    logging.info("Consulting id: %s/%s", i, m)
    try:
        with Entrez.esummary(db = "gds", id = id) as handle:
            metad = Entrez.read(handle)[0]
        return(metad)
    except:
        if not tries: logging.info("Error: %s", id); return({'Id':id})
        return(get_sum(id,i,m,tries-1))

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="Path to file with GSM list")
parser.add_argument("-o", "--out", help="File to print output")
args = parser.parse_args()

start = time()
with open(args.path) as file:
    gsms = [gsm.strip('\n') for gsm in file.readlines()]

m = len(gsms)
summaries = [get_sum(acc, i+1,m) for i,acc in enumerate(gsms)]

pd.DataFrame.from_dict(summaries).to_csv(args.out, sep='\t')
logging.info("Done in %s", format_time(time()-start))