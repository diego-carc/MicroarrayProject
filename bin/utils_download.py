import os
from ftplib import FTP, error_perm, error_temp
import logging
import GEOparse
import re
import psutil
import gzip
from shutil import copyfileobj

class FTPConnectionError(Exception):
    pass

class CheckSumError(Exception):
    pass

class DiskUsageError(Exception):
    pass

class File:
    
    def __init__(self, filePath):
        self.filePath = os.path.abspath(filePath)
        self.total_size = None
        
    def exists(self):
        return(os.path.isfile(self.filePath))
    
    def is_corrupted(self):
        if not self.exists() or (os.path.getsize(self.filePath) - self.total_size):
            return(True)
        return(False)
    
    def remove(self):
        if self.exists(): os.remove(self.filePath)
        else: logging.warning(f"{self.name} does not exist.")
    
    def is_gzip(self):
        if not self.exists(): 
            logging.warning(f"{self.filePath} does not exist")
            return 
        
        with open(self.filePath, "rb") as f:
            return True if f.read(2) == b'\x1f\x8b' else False
    
    def gunzip(self):
        if not self.is_gzip(): 
            logging.warning(f"{self.filePath} is not gzip")
            return 

        with gzip.open(self.filePath) as gz:
            with open(temp := self.filePath.strip(".gz"), "wb") as out:
                copyfileobj(gz, out)
        self.remove()
        self.filePath = temp


class GEOFile(File):
    host = "ftp.ncbi.nlm.nih.gov"
    
    def __init__(self, local_path, geo_path):
        self.geo_path = geo_path
        self.name = os.path.basename(self.geo_path)
        super().__init__(os.path.join(local_path, self.name))

    def is_downloaded(self):
        return super().exists()
    
    def is_corrupted(self):
        return super().is_corrupted()
    
    def remove(self):
        if self.is_downloaded(): os.remove(self.filePath)
        else: logging.warning(f"{self.name} has not been downloaded yet")

    def connect_to_ftp(self, tries=10):
        logging.info(f"Connecting to {self.host}...")
        try:
            ftp = FTP(self.host)
            ftp.login()
            return(ftp)
        except:
            if tries: return(self.connect_to_ftp(tries=tries-1))
            logging.warning(f"Unnable to connect to {self.host}")
            raise FTPConnectionError
        
    def __download__(self, ftp:FTP):
        # Try to download
        self.total_size = ftp.size(self.geo_path)
        if self.total_size >= psutil.disk_usage("/").free: raise DiskUsageError
        logging.info(f"{self.name} total size: {self.total_size}")
        logging.info(f"Downloading {self.name} to {self.filePath}...")
        with open(self.filePath, "wb") as f:
            ftp.retrbinary(f"RETR {self.geo_path}", f.write)
        
    def download(self, ftp:FTP=None, close:bool=True, tries:int=10):
        if self.is_downloaded(): 
            logging.info(f"{self.name} has been succesfully downloaded")

            if close and ftp: ftp.quit(); ftp = None
            return(ftp)

        try:
            # Try to connect
            if not ftp: ftp = self.connect_to_ftp()
            
            # Try to download
            self.__download__(ftp)

            # Check integrity
            logging.info("Checking integrity...")
            if not self.total_size: self.total_size = 0
            if self.is_corrupted(): raise CheckSumError

        except error_perm:
            logging.warning(f"File {self.geo_path} not found in {self.host}")
            tries = 0
        except FTPConnectionError: 
            logging.warning(f"Unable to download {self.name}") 
            ftp = None
        except CheckSumError: 
            logging.warning(f"{self.name} file size does not match expected size from {self.geo_path}")                
            self.remove()
        except EOFError: 
            logging.warning(f"FTP connection failed attempting to download {self.name}")
            ftp = None
            self.remove()
        except error_temp as e:
            logging.warning(f"{e} ocurred attempting to download {self.name}")
            ftp = None
            self.remove()
        except DiskUsageError:
            logging.critical(f"No available storage in disk to download {self.name}. Aborting.")
            exit()
        if tries: return(self.download(ftp, close, tries-1))
        logging.warning(f"The following file wont be downloaded: {self.name}")
        if close and ftp: ftp.quit(); ftp = None
        return(ftp)

def is_processable(row):
    """
    
    """
    smfile = row.GEO_file_path
    man = row.GPL_manufacturer
    dist = row.GPL_distribution
    procesable = False
    if type(smfile) == float or type(man) == float or type(dist) == float: pass
    elif re.search(r'norm_RMA', smfile, flags=re.IGNORECASE): pass
    elif dist == "commercial" or dist == "custom-commercial": 
        if man == "Affymetrix" and re.search(r"\.cel|\.exp", smfile, flags=re.IGNORECASE) : procesable = True
        if man == "Agilent" and re.search(r"\.txt", smfile, flags=re.IGNORECASE): procesable = True
        if man == "NimbleGen" and re.search(r"\.pair|\.xys", smfile, flags=re.IGNORECASE): procesable = True
        if man == "Illumina" and re.search(r"\.idat", smfile, flags=re.IGNORECASE): procesable = True
    return(procesable)

def get_GEO(geo_id:str, path:str, remove:bool=False, tries:int=10):
    '''
    Calls GEOparse to parse GSE SOFT file and retrieve a GSE object

    Args:
        geo_id(str) : GEO accession id
        path(str): Path to dir to download SOFT file
        remove(bool): Wether to delete downloaded SOFT file
        tries(int): Number of tries before returning None
    Return:
        GSE object or None
    '''
    geo_type = {"GPL":"platforms", "GSM":"samples", "GSE":"series"}
    geo_path = f'geo/{geo_type[geo_id[:3]]}/{geo_id[:-3]}nnn/{geo_id}/soft/{geo_id}_family.soft.gz'
    geo_file = GEOFile(path, geo_path)

    try:
        try:
            logging.info(f"Attempting to parse {geo_id}...") 
            geo = GEOparse.get_GEO(geo=geo_id, destdir=path, silent=True)
            if geo_id[:3] == "GSM" or geo_id[:3] == "GPL": geo_file.filePath = os.path.join(path, f"{geo_id}.txt")
            logging.info(f"{geo_id} GEOparse succesful")
        except:
            logging.warning("An error ocurred with GEOparse. Attempting to download via ftp...")
            geo_file.download()
            geo = GEOparse.get_GEO(filepath=geo_file.filePath, silent=True)
            logging.info(f"{geo_id } GEOparse succesful")
        if remove: geo_file.remove()
    except:
        if tries: return(get_GEO(geo_id, path, remove, tries-1))
        logging.warning(f"Unable to parse {geo_id}") ; return(None)
    return(geo)

def save_table(geo:GEOparse.GPL|GEOparse.GSE|GEOparse.GSM, path:str, **kwargs):
    try:
        geo.table.to_csv(f"{path}", **kwargs)
        logging.info(f"{geo.name} table saved")
    except:
        logging.warning(f"Unable to save {geo} table")
    

