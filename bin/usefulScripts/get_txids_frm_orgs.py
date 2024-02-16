from Bio import Entrez
import pandas as pd
from time import time
Entrez.email = "diegocar@lcg.unam.mx"
Entrez.api_key=  "34677fdcfd2f0659a7f9ee05ab6e44704f09"
from functions import format_time
import argparse
import logging
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)

def get_org_id(org, i, m,tries=10):
    logging.info("Consulting org %s/%s", i, m)
    try:
        with Entrez.esearch(db="taxonomy", term=org) as handle:
            record = Entrez.read(handle)
            return(record['IdList'][0])
    except:
        if tries: return(get_org_id(org, i, m, tries -1))
        logging.info("Error: %s", org)
        return("None")
     


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="Path to file with Orgs list")
parser.add_argument("-o", "--out", help="Path + file name to print output")
args = parser.parse_args()

start = time()
with open(args.path) as file:
    orgs = [org.strip('\n') for org in file.readlines()]

m = len(orgs)
orgs_ids = {org : get_org_id(org, i+1,m) for i,org in enumerate(orgs)}
pd.Series(orgs_ids).to_csv(args.out, header=False, sep="\t")
logging.info("Done in %s", format_time(time()-start))