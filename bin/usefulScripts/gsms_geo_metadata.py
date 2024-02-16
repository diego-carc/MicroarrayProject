import GEOparse as g
import argparse
import logging
import os
from multiprocessing import Pool
import pandas as pd

def get_gsm_object(gsm_id, path='./'):
  try:
    try: gsm = g.get_GEO(geo=gsm_id, destdir=path, silent=True)
    except: 
       link = f'"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={gsm_id}&form=text&view=full"'
       os.system(f"wget {link} -O {gsm_id}.txt") 
       gsm = g.get_GEO(filepath=f"./{gsm_id}.txt", silent=True)
    os.remove(f'{path}/{gsm_id}.txt')
    return(gsm)
  except:
    logging.warning("No se pudo parsear el GSM %s", gsm_id)
    return(None)
  
def get_gsm_metadata(gsm_id):
   logging.info("Parsing %s", gsm_id)
   if not (gsm := get_gsm_object(gsm_id)): return(gsm)
   return({k: ';'.join(v) for k,v in  gsm.metadata.items()})

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gsms", help="Path to file with GSMs list")
parser.add_argument("-o", "--out", help="Path to outputfile")
parser.add_argument("-c", "--cores", help="Cores", default=1, type=int)
args = parser.parse_args()

with open(args.gsms) as file:
    gsms  = list({gsm.strip('\n') for gsm in file.readlines()})

with Pool(args.cores) as p:
  metadata_sync = p.map(get_gsm_metadata, gsms)

pd.DataFrame.from_dict([i for  i in metadata_sync if i]).to_csv(args.out,sep = '\t', header=True, index=False)