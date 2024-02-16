"""

"""
import os
import pandas as pd
import numpy as np
import re

class ColumnNameError(Exception):
    pass

class EncodingError(Exception):
    pass
     
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
    pass

class Agilent(Parser):
    
    def __init__(self, filePath):
        for p,d in [("ProbeName", "gMedianSignal"), ("ID", "Raw intensity (med) {532}"), ("Name", "F532 Median"),
                    ("Name", "Signal Median"), ("ID", "F635 Median")]:
            try: 
                super().__init__(filePath, p, d)
            except ColumnNameError: pass 

    def parseFile(self):
       return super().parseFile()
    
class NimbleGen(Parser):
    
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

def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values, axis=0), index=df.index, columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)