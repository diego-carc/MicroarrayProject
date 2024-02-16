"""
Select channel
"""

# Import modules
from os.path import abspath, join, basename
from os import makedirs
from utils_normalization import *
from utils_download import File
import pandas as pd
from time import time
from utils_query import format_time
import argparse
import logging
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)

# Parse args
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Path to file downloadPaths.tsv")
parser.add_argument("-o", "--output", help="Path to directory to save normalized files")
parser.add_argument("-r", "--raw", help="Use this flag to save raw files", action="store_true")
args = parser.parse_args()
log = f"""Running normalization with the following parameters:
\tinput: {args.input}
\toutput: {args.output}
\traw: {args.raw}
"""
logging.info(log)

# Function
def handle(filePath:str, GPL_manufacturer:str, message:str=None):
    logging.info(message)
    manufacturers = {"Agilent": Agilent, "NimbleGen":NimbleGen}
    try: data = manufacturers.get(GPL_manufacturer)(filePath).data
    except Exception as error: logging.info(f"Error parsing {filePath}: {error}"); data = pd.DataFrame.from_dict({})
    return data

# Read input
logging.info(f"Reading {args.input}...")
data = pd.read_csv(args.input, sep='\t')
logging.info("Done!")
start = time()

# Select processable and downloaded
logging.info("Filtering...")
data = data[data.is_processable & data.is_downloaded]
logging.info("Done!")

# Get destiny paths
logging.info("Parsing local paths...")
dirs =  ["Specie", "TaxId", "GPL_geo_accession"] 
getLocalPath = lambda row : abspath(join(args.output, *[dest if d != "TaxId" else f"TaxId_{dest}" 
                                           for d in dirs if (dest := str(row.get(d)).replace(" ", '_')) != "None"]))
data["Normalized_path"] = data.apply(getLocalPath, axis=1)
logging.info("Done!")

# Get file names
logging.info("Parsing file names...")
data["Normalized_name"] = data.GPL_geo_accession.apply(lambda GPL: f"{GPL}_normalized.tsv")
if args.raw: data["Raw_name"] = data.GPL_geo_accession.apply(lambda GPL: f"{GPL}_raw.tsv")
logging.info("Done!")

# Gunzip files
data["Files"] = data.apply(lambda row : File(join(row.Local_path, row.File_name)), axis=1)
data.Files.apply(lambda f: f.gunzip())
data.File_name = data.Files.apply(lambda f: basename(f.filePath))


# Normalize
logging.info("Starting normalization...")
totalDirectories = data.Local_path.unique().shape[0]
for i,(path, tab) in enumerate(data.groupby("Local_path")):
    
    logging.info(f"Parsing files in {path}")
    logging.info(f"Directory {i+1}/{totalDirectories}")
    parsedFiles = [parsed for i,(indx,row) in enumerate(tab.iterrows()) 
                   if not (parsed:=handle(row.Files.filePath, row.GPL_manufacturer,
                                      f"Parsing file: {i+1}/{tab.File_name.shape[0]}")).empty]
    
    logging.info(f"Merging data in {path}...")

    if not parsedFiles:
        logging.warning(f"Could not parse files in: {path}")
        continue

    merged = pd.concat(parsedFiles)
    if args.raw: 
        makedirs(tab.iloc[0,:].Normalized_path, exist_ok=True)
        logging.info(f"Saving raw data...")
        merged.to_csv(join(tab.iloc[0,:].Normalized_path, tab.iloc[0,:].Raw_name), sep='\t', index=False)

    logging.info("Applying quantile normalization...")
    try:
        makedirs(tab.iloc[0,:].Normalized_path, exist_ok=True)
        norm = quantile_normalize(merged)
        logging.info(f"Saving normalized data...")
        norm.to_csv(join(tab.iloc[0,:].Normalized_path, tab.iloc[0,:].Normalized_name), sep='\t', index=False)  
    except: 
        logging.warning(f"Could not normalize files in {path}")  

# Report
logging.info("Saving normalization paths...")
data = data.drop("Files", axis=1)
data.to_csv(join(args.output, f"normalizationPaths.tsv"),sep='\t', index=False)
logging.info(f"Normalization terminated in {format_time(time() - start)}")