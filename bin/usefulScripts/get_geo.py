import GEOparse as g
import os
import argparse
import pandas as pd

def get_geo(geo, i, m, tries = 10):
    print(f"{i}/{m}") 
    try:
        link = f'"http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={geo}&form=text&view=full"'
        os.system(f"wget {link} -O {geo}.txt") 
        o = g.get_GEO(filepath=f"./{geo}.txt", silent=True)
        os.remove(f"./{geo}.txt")
        return(o)

    except: 
        print(f"Error in {geo}")
        if tries: return(get_geo(geo, i, tries-1))
        return(None)

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--geos", help="Path to file with GEO accessions")
parser.add_argument("-o", "--out", help="Path to out file")
args = parser.parse_args()

with open(args.geos) as file:
    geo_acc = list({ga.strip('\n') for ga in file.readlines()})

m = len(geo_acc)

# Get geo objects
geo_objc = {gpl: o for i,gpl in enumerate(geo_acc) if (o:=get_geo(gpl, i, m))}

# Parse metadata
geo_meta = {gpl : {i:';'.join(j) for i,j in o.metadata.items()}  for gpl,o in geo_objc.items()}

pd.DataFrame.from_dict(geo_meta).T.reset_index().rename(columns={"index":"GPL"}).to_csv(f"{args.out}", sep='\t', index=False)