"""

"""
import pandas as pd
import numpy as np
import re
import logging
from sys import stdout
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='[%m/%d/%Y %I:%M:%S %p]', level=logging.INFO, stream=stdout)

class ColumnNameError(Exception):
    pass

class EncodingError(Exception):
    pass

class MicroarrayFileFormatError(Exception):
    pass

class Microarray():
    import pandas as pd
    from utils_download import File
    import re
    import os
    
    
    def __init__(self, dataFilePath:str, annotFilePath:str, GSM:str=None,
                 oneChannel:bool=True, annotCols:set={}, dataCols:set={}):
        
        # Sample identifier
        self.GSM =  GSM if GSM else self.re.search(r"(GSM\d+)", dataFilePath).group()
        logging.debug(f"Microarray sample name is: {self.GSM}")

        # Channel count
        self.oneChannel = oneChannel

        # File with expression data
        self.dataFile = self.File(dataFilePath)
        logging.debug(f"Microarray data file is: {self.dataFile.filePath}")

        # File with gene annotation data
        self.annotFile = self.File(annotFilePath)
        logging.debug(f"Microarray annot file is: {self.annotFile.filePath}")

        # Columns with expression and gene data
        self.annotCols = annotCols
        self.dataCols = dataCols
        logging.debug(f"Columns are {self.annotCols} and {self.dataCols}")

        # DataFrames with Microarray data
        self.expressionData = self.parseRawData()
        logging.debug(f"Raw data parsed")
        self.annotationData = None
        logging.debug(f"Annot data parsed")
    
        logging.debug("Microarray init succesful")
    
    def parseRawData(self, **kwargs) -> pd.DataFrame:
        """ Parse a tab delimited file

        _extended summary_[#_unique ID_]_


        Returns:
            pd.DataFrame: DataFrame containing the parsed file

        References:
            .. [#_unique ID_] _pubmed abbr journal title_ _vol_:_page or e-article id_ (_year_) https://doi.org/_doi_
            .. [#_unique ID_] _first-author first-name last-name_ _book title_ (_year_) ISBN:_ISBN_ _http link_
            .. [#_unique ID_] _article title_ _conference_ (_year_) _http link_"""
        try:
            file = self.pd.read_csv(self.dataFile.filePath, sep='\t', **kwargs)
        except UnicodeDecodeError:
            file = self.pd.read_csv(self.dataFile.filePath, sep='\t', encoding="latin-1", **kwargs)
        finally:
            logging.debug(f"{self.dataFile.filePath} parsed")
            return file

class Agilent(Microarray):
    def __init__(self, dataFilePath:str, annotFilePath:str="", oneChannel:bool=True):
        """ Initiates an Agilent Microarray object

        _extended summary_[#_unique ID_]_

        
        Arguments:
            dataFilePath (str): _description_
            annotFilePath (str): _description_. Defaults to "".
            oneChannel (bool): _description_. Defaults to True.

        References:
            .. [#_unique ID_] _pubmed abbr journal title_ _vol_:_page or e-article id_ (_year_) https://doi.org/_doi_
            .. [#_unique ID_] _first-author first-name last-name_ _book title_ (_year_) ISBN:_ISBN_ _http link_
            .. [#_unique ID_] _article title_ _conference_ (_year_) _http link_"""
        
        super().__init__(dataFilePath = dataFilePath,
                         annotFilePath = annotFilePath,
                         oneChannel = oneChannel, 
                         annotCols={"SystematicName", "GeneName"},
                         dataCols={"rMedianSignal", "gMedianSignal"})
        logging.debug("Agilent Microarray init succesful")


    def parseRawData(self) -> pd.DataFrame:
        """Reads Agilent Microarray supplementary files

        _extended summary_[#_unique ID_]_


        Raises:
            MicroarrayFileFormatError: _description_

        Returns:
            pd.DataFrame: _description_

        References:
            .. [#_unique ID_] _pubmed abbr journal title_ _vol_:_page or e-article id_ (_year_) https://doi.org/_doi_
            .. [#_unique ID_] _first-author first-name last-name_ _book title_ (_year_) ISBN:_ISBN_ _http link_
            .. [#_unique ID_] _article title_ _conference_ (_year_) _http link_"""
        
        # Check if file has Agilent format
        if self.dataFile.readChars(3) != "TYP": 
            raise MicroarrayFileFormatError(f"{self.dataFile.filePath} might not be an Agilent File")
        
        # Parse file
        rawData = super().parseRawData(skiprows=range(9)).rename(columns=lambda col : col.replace(' ',''))
        logging.debug(f"Raw data parsed")

        # Filter control probes & select columns
        rawData = rawData[rawData.ControlType == 0][list(set(rawData.columns) & (self.annotCols | self.dataCols))]
        logging.debug(f"Control probes filtered and colums selected")
        
        # Select annotation column
        if not (annot:="SystematicName") in rawData.columns: annot = "GeneName"
        elif "GeneName" in rawData.columns: rawData = rawData.drop("GeneName", axis=1) 
        
        # Mean redundancies and rename (only two channel)
        rawData = rawData.groupby(annot).mean().rename(columns=lambda col : self.GSM if (self.oneChannel) else  f"{self.GSM}_{col[0]}") 
        logging.debug("Redundancies meaned")
        return rawData
    
    def parseAnnotData(self):
        pass

class Affymetrix(Microarray):
    def __init__(self, dataFilePath: str, annotFilePath: str, GSM: str = None, oneChannel: bool = True, annotCols: set = {}, dataCols: set = {}):
        super().__init__(dataFilePath, annotFilePath, GSM, oneChannel, annotCols, dataCols)
    
class NimbleGen(Microarray):
    def __init__(self, dataFilePath: str, annotFilePath: str='', oneChannel: bool = True):
        # Determine file format
        self.frmt = self.os.path.splitext(dataFilePath)[-1].lower()
        if self.frmt != ".pair" and self.frmt != ".xys": 
            raise MicroarrayFileFormatError(f"{dataFilePath} migth not a NimbleGen file")
        logging.debug(f"{dataFilePath} has {self.frmt} format")

        # Choose proper columns
        if self.frmt == ".pair":
            annotCols = {"SEQ_ID"}
            dataCols= {"PM"}
        elif self.frmt == ".xys":
            annotCols = {"X", "Y"}
            dataCols = {"SIGNAL"}
        logging.debug(f"{dataFilePath} columns are {annotCols | dataCols}")

        # Init Microarray
        super().__init__(dataFilePath = dataFilePath,
                         annotFilePath = annotFilePath,
                         oneChannel = oneChannel,
                         annotCols = annotCols,
                         dataCols = dataCols)
        logging.debug("NimbleGen init succesful")

    def parseRawData(self) -> pd.DataFrame:
        # Count rows to skip
        commentLines = 1 if self.dataFile.readChars(1) == "#" else 0
        logging.debug(f"{self.dataFile.filePath} has {commentLines} comment lines")
        
        # Read file
        rawData = super().parseRawData(skiprows = range(commentLines))
        logging.debug(f"{self.dataFile.filePath} read succesfully")

        # Delete RANDOM probes
        if self.frmt == ".pair":
            rawData = rawData[rawData.GENE_EXPR_OPTION != "RANDOM"]
            logging.debug("RANDOM probes deleted")
        else: pass

        # Select Columns
        rawData = rawData[list(set(rawData.columns) & (self.annotCols | self.dataCols))]
        logging.debug(f"Columns are {rawData.columns}")
        if self.frmt == ".pair":
            annot = list(self.annotCols)[0]
        elif self.frmt == ".xys":
            annot = "COORDS"
            rawData["COORDS"] = rawData.apply(lambda row: '_'.join([str(int(row.X)), str(int(row.Y))]),axis=1)
            rawData = rawData.drop(["X", "Y"], axis=1).dropna()
        logging.debug(f"Annot column is {annot}")

        # Mean redundancies and rename
        return rawData.groupby(annot).mean().rename(columns=lambda col: self.GSM)

class GenePix(Microarray):
    def __init__(self, dataFilePath, annotFilePath):
        super().__init__(dataFilePath, annotFilePath)

class Parser:

    def __init__(self, filePath, probeCol, dataCol):
        self.filePath = filePath
        self.probeCol = probeCol
        self.dataCol = dataCol
        self.start = self.findLine() 
        if self.start == None: raise ColumnNameError(f"{self.probeCol} and {self.dataCol} failed as column names")
        self.GSM = re.search(r"(GSM\d+)", filePath).group()
        self.data = self.parseFile()
    
    def parseFile(self):
        """
        Reads a tab separated file into a pandas DataFrame
        """
        try: 
            data = pd.read_csv(self.filePath, sep='\t', skiprows=range(self.start), usecols=[self.probeCol, self.dataCol])
        except: 
            data = pd.read_csv(self.filePath, sep='\t', skiprows=range(self.start), usecols=[self.probeCol, self.dataCol], encoding="latin-1")

        if data.dtypes[self.dataCol] == "O": data[self.dataCol] = data[self.dataCol].apply(lambda strVal: float(strVal.replace(',', '.')))

        return self.meanRedundancies(data) 

    def meanRedundancies(self, dataFrame):
        """
        Groups the data by the repeated ProbeNames and returns the mean

        """
        grouped = {pName : data[self.dataCol].mean() for pName, data in dataFrame.groupby(self.probeCol)}
        grouped = pd.DataFrame.from_dict({"ProbeName":grouped.keys(), self.GSM:grouped.values()})
        return grouped.set_index("ProbeName")
    
    def findLine(self, encoding="utf-8"):
        """
        Iterates over the lines of a tab separated file to find the header line
        """
        try:
            with open(self.filePath, encoding=encoding) as file:
                for i,line in enumerate(file.readlines()):
                    if self.dataCol in line and self.probeCol in line:
                        return i
                return None
        except: 
            if encoding == "latin-1": raise EncodingError(f"Could not open {self.filePath}")
            return(self.findLine("latin-1"))

    
class Affymetrix(Parser):
    
    def __init__(self, filePath, probeCol, dataCol):
        super().__init__(filePath, probeCol, dataCol)

    
class oldNimbleGen(Parser):
    
    def __init__(self, filePath, probeCol = "PROBE_ID", dataCol="PM"):
        super().__init__(filePath, probeCol, dataCol)
    
    def parseFile(self):
        return super().parseFile()
    
class GenePix(Parser):

    def __init__(self, filePath, probeCol, dataCol):
        super().__init__(filePath, probeCol, dataCol)

class Illumina(Parser):
    pass

class Miscelaneous(Parser):
    pass

class MicroarrayBatch():
    
    def quantile_normalize(self):
        df = self
        """
        input: dataframe with numerical columns
        output: dataframe with quantile normalized values
        """
        df_sorted = pd.DataFrame(np.sort(df.values, axis=0), index=df.index, columns=df.columns)
        df_mean = df_sorted.mean(axis=1)
        df_mean.index = np.arange(1, len(df_mean) + 1)
        df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
        return(df_qn)