"""
"""
import pandas as pd
import numpy as np
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
    
    
    def __init__(self, dataFilePath:str, annotFilePath:str, GSM:str='',
                 oneChannel:bool=True, annotCols:set=set(), dataCols:set=set()):
        
        # Sample identifier
        if not GSM and not (gsm := self.re.search(r"(GSM\d+)", dataFilePath)):
            raise MicroarrayFileFormatError(r"Please provide a GSM to identify the sample or\
                                            supply a filename matching GSM\d+ pattern")
        self.GSM =  GSM if GSM else gsm.group()
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
    
    def parseAnnotData(self):
        pass

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
    def __init__(self, dataFilePath: str, annotFilePath: str, GSM: str = '', oneChannel: bool = True, annotCols: set = {}, dataCols: set = {}):
        super().__init__(dataFilePath, annotFilePath, GSM, oneChannel, annotCols, dataCols)
    
class NimbleGen(Microarray):
    def __init__(self, dataFilePath: str, annotFilePath: str='', oneChannel: bool = True):
        # Determine file format
        self.frmt = super().os.path.splitext(dataFilePath)[-1].lower()
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

class GenPix(Microarray):
    from pandas.core.api import DataFrame as DataFrame

    def __init__(self, dataFilePath: str, annotFilePath: str = '', GSM: str = '', oneChannel: bool = True):
        """Initiate a GenPix Microarray

        _extended summary_[#_unique ID_]_

        .. math:: _LaTeX formula_

        Arguments:
            dataFilePath (str): _description_
            annotFilePath (str): _description_
            GSM (str): _description_. Defaults to ''.
            oneChannel (bool): _description_. Defaults to True.
            annotCols (set): _description_. Defaults to set().
            dataCols (set): _description_. Defaults to set().

        References:
            .. [#_unique ID_] _pubmed abbr journal title_ _vol_:_page or e-article id_ (_year_) https://doi.org/_doi_
            .. [#_unique ID_] _first-author first-name last-name_ _book title_ (_year_) ISBN:_ISBN_ _http link_
            .. [#_unique ID_] _article title_ _conference_ (_year_) _http link_"""
        
        # Call superclass init
        logging.debug("Initiation GenPix")
        super().__init__(dataFilePath, annotFilePath, GSM, oneChannel,
                          annotCols = {"Name", '"Name"', "ID", '"ID"'},
                        dataCols = {"F532 Median", "F635 Median", '"F532 Median"', '"F635 Median"'})
        logging.debug("GenPix Initiation succesflu!")
    
    def parseRawData(self) -> DataFrame:
        """ Parse raw gpr data files

        _extended summary_[#_unique ID_]_

        .. math:: _LaTeX formula_

        Raises:
            MicroarrayFileFormatError: _description_

        Returns:
            DataFrame: _description_

        References:
            .. [#_unique ID_] _pubmed abbr journal title_ _vol_:_page or e-article id_ (_year_) https://doi.org/_doi_
            .. [#_unique ID_] _first-author first-name last-name_ _book title_ (_year_) ISBN:_ISBN_ _http link_
            .. [#_unique ID_] _article title_ _conference_ (_year_) _http link_"""

        # Check if file has GPR format
        if self.dataFile.readChars(3) != "ATF":
            raise MicroarrayFileFormatError(f"{self.dataFile.filePath} might not be a GenPix File")
        
        # Count comment lines
        logging.debug("Counting comment lines")
        commLines = self.countCommentLines()
        logging.debug(f"{self.GSM} has {commLines} comment lines")

        # Parse table
        logging.debug(f"Parsing {self.GSM}")
        rawData =  super().parseRawData(skiprows=range(commLines))
        logging.debug(f"File parsed. Dim: {rawData.shape}")

        # Rename dyes
        rawData = rawData.rename(columns=lambda col: col.replace("647", "635").replace("555", "532"))
        logging.debug("Colums renamed to use F635 and F532 nomenclature")

        # Select columns
        rawData = rawData[list(set(rawData.columns) & (self.annotCols | self.dataCols))]
        logging.debug(f"Columns selected: {rawData.columns}")

        # Select annotation column
        ## Use column with the less number of NaN
        annotColsPresent = self.annotCols & set(rawData.columns)
        annot =  sorted(annotColsPresent, key=lambda col: sum(rawData[col].isna()))[0]
        rawData = rawData.drop(annotColsPresent - {annot},axis=1)
        logging.debug(f"Annotation col selected: {annot}")

        # Mean and Rename
        chName = lambda col : f"{self.GSM}_g" if "532" in col else f"{self.GSM}_r"
        rawData = rawData.groupby(annot).mean().rename(columns=lambda col : self.GSM if (self.oneChannel) else chName(col))
        logging.debug("Redundancies meaned and columns renamed")

        return rawData
    
    def countCommentLines(self) -> int:
        """ Count comment lines to skip

        _extended summary_[#_unique ID_]_

        .. math:: _LaTeX formula_

        Returns:
            int: _description_

        References:
            .. [#_unique ID_] _pubmed abbr journal title_ _vol_:_page or e-article id_ (_year_) https://doi.org/_doi_
            .. [#_unique ID_] _first-author first-name last-name_ _book title_ (_year_) ISBN:_ISBN_ _http link_
            .. [#_unique ID_] _article title_ _conference_ (_year_) _http link_"""
        try: 
            with open(self.dataFile.filePath) as file:
                temp = file.readlines()
        except UnicodeDecodeError: 
            with open(self.dataFile.filePath, encoding="latin-1") as file:
                temp = file.readlines()
        finally:
            # Check if line has intersection with annotation columns
            for i, line in enumerate(temp):
                if self.annotCols & set(line.split("\t")): return i
            return -1


class MicroarrayBatch():
    import numpy as np

    
    def quantile_normalize(self, df):
        """
        input: dataframe with numerical columns
        output: dataframe with quantile normalized values
        """
        df_sorted = pd.DataFrame(np.sort(df.values, axis=0), index=df.index, columns=df.columns)
        df_mean = df_sorted.mean(axis=1)
        df_mean.index = self.np.arange(1, len(df_mean) + 1)
        df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
        return(df_qn)