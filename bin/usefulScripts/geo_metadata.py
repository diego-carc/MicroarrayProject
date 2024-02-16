"""
Calls GEOparse to build a DataFrame containing the metadata for a list of GPLs accession numbers
"""
from GEOparse import get_GEO
from functions import format_time
from time import time
import logging
from os import remove
import argparse
import pandas as pd
 
def get_geo_object(ID, path, i,m, tries = 10):
  logging.info("Parsing GPL %s/%s",i,m)
  try:
    geo = get_GEO(geo=ID, destdir=path, silent=True)
    remove(f'{path}/{ID}.txt')
    return(geo)
  except:
    if not tries: return(None); logging.info("Error: %s", ID)
    return(get_geo_object(ID, path,i,m, tries-1))
  
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)

# Parse args
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="Path to file with GEO accession numbers")
parser.add_argument("-o", "--out", help="Output file path name")
args = parser.parse_args()

# Read GPL list
with open(args.path) as file:
    geos = [i.strip('\n') for i in file.readlines()]

# For each id in GPL list get gpl object
start = time()
m = len(geos)
a = {gpl : get_geo_object(gpl,'./', i+1, m) for i,gpl in enumerate(geos)}
b = {k : {i: ';'.join(j) for i,j in v.metadata.items()} if v else {} for k,v in a.items()}
pd.DataFrame.from_dict(b).T.reset_index().rename(columns={'index':'GPL'}).to_csv(args.out, sep='\t')

logging.info("Done in %s", format_time(time() - start))