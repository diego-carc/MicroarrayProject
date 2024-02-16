from utils_download import GEOFile
import pandas as pd 
import os


out = "/export/storage/users/diegocar/abasySpecieFiles/GPL_supplementary_files"

data = pd.read_csv("/export/storage/users/diegocar/Test/Query/2312/07_getGeoTest/getGeoTest.tsv", sep='\t')

files = [GEOFile(out, path) for paths in data.GPL_supplementary_file if type(paths) != float for path in paths.split(";")]

total = len(files)
ftp = None
for i,file in enumerate(files):
    print(f"Downloading: {i+1}/{total}")
    try: 
        ftp = file.download(ftp = ftp, close=False)
    except: 
        ftp = None
        print(f"Error: {file.name}")
if ftp: ftp.close()

print("Done!")


