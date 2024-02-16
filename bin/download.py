'''
Download supplementary files
Download GEO data https://www.ncbi.nlm.nih.gov/geo/info/download.html
FTP README https://ftp.ncbi.nlm.nih.gov/geo/README.txt


'''
# Import modules
from os import listdir, makedirs
from os.path import join, abspath, split
import argparse
import pandas as pd
from utils_download import *
from utils_query import format_time
import logging
from time import time
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO)

# Parse args
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Metadata table from query code")
parser.add_argument("-o", "--outdir", help=f"Directory to download files. Def: {(temp:='./')}", default=temp)
parser.add_argument("-m", "--misc", help="Wether to download miscelaneous (no processable) files", action="store_true")
parser.add_argument("-t", "--taxrel", help="Path to taxonomy relation dataframe to group by specie", default=None)
args = parser.parse_args()
log = f"""Running download with the following parameters:
\tinput: {args.input}
\toutdir: {args.outdir}
\tmisc: {args.misc}
\ttaxrel: {args.taxrel}"""
logging.info(log)

# Track total time
total_time = time()

# Read table
logging.info(f"Reading data from {args.input}...")
data = pd.read_csv(args.input, sep='\t', usecols=["TaxId", "GPL_geo_accession", "GPL_manufacturer" , "GSM_supplementary_file",
                                                   "GSM_geo_accession", "GPL_distribution"])

# Merge taxonomy relations
if args.taxrel:
    logging.info("Merging taxonomy relations")
    data = data.merge(pd.read_csv(args.taxrel, sep='\t', usecols=["TaxId", "Specie"]))

# Make a row for each supplementary file
logging.info("Solving multiple supplementary files...")
files = [{**row, "GEO_file_path" : f} for i,row in data.iterrows() for f in row.GSM_supplementary_file.split(';')]
files = [{k:v for k,v in file.items() if k != "GSM_supplementary_file"} for file in files]
data = pd.DataFrame.from_dict(files)
logging.info("Done!")

# Determine processability
logging.info("Classifying supplementary files...")
data["is_processable"] = data.apply(is_processable, axis=1)
logging.info("Done!")

# Get local paths
logging.info("Parsing local paths...")
dirs =  ["Specie", "TaxId", "GPL_geo_accession"] 
getLocalPath = lambda row : abspath(join(args.outdir if row.is_processable else join(args.outdir, "Miscelaneous"), 
                                         *[dest if d != "TaxId" else f"TaxId_{dest}" 
                                           for d in dirs if (dest := str(row.get(d)).replace(" ", '_')) != "None"]))

data["Local_path"] = data.apply(getLocalPath, axis=1)
logging.info("Done!")

# Create file objects
logging.info("Preparing download...")
processable = data[(data.is_processable | args.misc)]
processable["Files"] = processable.apply(lambda row: GEOFile(row.Local_path, row.GEO_file_path), axis=1)
logging.info("Done!")

# Download supplementary files
supp_time = time()
ftp = None
total = processable.shape[0]
for i,f in enumerate(processable[processable.GEO_file_path != "NONE"].Files):
    logging.info(f"Downloading file {i+1}/{total}..")
    makedirs(split(f.filePath)[0], exist_ok=True)
    ftp = f.download(ftp, close=False)
if ftp: ftp.quit()
logging.info("Download completed!")

# Count downloaded
processable["File_name"] = processable.Files.apply(lambda f: f.name)
processable["is_downloaded"] = processable.Files.apply(lambda f: f.is_downloaded())

# Get GPL tables
logging.info("Downloading GPL tables...")
gpl_time = time()
gpl_dir = join(args.outdir, "GPL_Tables")
makedirs(gpl_dir, exist_ok=True)
for gpl in data.GPL_geo_accession.unique():
    gpl_table = GEOFile(gpl_dir, join("temp", f"{gpl}.txt"))
    if gpl_table.is_downloaded(): continue
    save_table(get_GEO(geo_id=gpl, path=gpl_dir, remove=True), gpl_table.filePath, sep='\t', index=False)
gpl_end = time()
logging.info("GPL tables download!")

# Report
logging.info(f"{len(listdir(gpl_dir))} GPL tables downloaded in {format_time(gpl_end - gpl_time)}")
logging.info(f"{data[data.is_processable].shape[0]} files were classified as processable")
logging.info(f"{data[~data.is_processable].shape[0]} files were classified as miscelaneous (no processable)")
logging.info(f"{sum(processable.is_downloaded)}/{total} files were downloaded in {format_time(time()-supp_time)}!")
logging.info(f"Recovery: {sum(processable.is_downloaded)*100/total}%")
logging.info(f"Total time: {format_time(time()-total_time)}!")
logging.info(f"Saving Paths table to {(temp :=  join(args.outdir, "downloadPaths.tsv"))}")
processable = processable.drop(["Files", "GEO_file_path"], axis=1)
processable.to_csv(temp, index=False, sep='\t')
logging.info("Done!")