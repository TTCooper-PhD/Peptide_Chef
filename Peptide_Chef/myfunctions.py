### Data Handling
import numpy as np
import pandas as pd
from pandas import ExcelWriter
import openpyxl
import re
import pyteomics
from pyteomics import fasta, parser, electrochem, mass
from itertools import combinations
from urllib.request import urlretrieve
import gzip
#dfdfdf
#Statistics
import statsmodels.api as sm
from statsmodels.formula.api import ols

#Figure Generation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import offsetbox
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import NullFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import axes3d
from IPython.display import Image, display
import seaborn as sns
from adjustText import adjust_text
import glob
import bioinfokit
from bioinfokit import analys, visuz

#Venn Diagrams
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
%matplotlib inline

cmap = 'PRGn'
fmt='eps'
dpi=600

# Use if Fasta file is Local
def Cook_Book(Species=None,homebrew=True, takeout=False,url=None ):
    ingredients=list()
    if takeout is True:
        homebrew=False
        print('Downloading the FASTA file from url...')
        urlretrieve(url,'temp.fasta.gz')
        print('Unzipping...')
        with gzip.open('temp.fasta.gz', mode='rt') as gzfile:
            for info, contents in fasta.FASTA(gzfile):
                taste=list((info,contents))
                ingredients.append(taste)
                recipie=pd.DataFrame(ingredients,columns=['ID','Peptide'])
        print("Takeout is Done!")
    if homebrew is True:
        print("Downloading the FASTA file from local flle...")
        book = str(Species+".fasta")
        recipie=pd.DataFrame()
        print("Serving up a homebrew...")
        with fasta.read(book) as menu:
            for info, contents in menu:
                taste=list((info,contents))
                ingredients.append(taste)
                recipie=pd.DataFrame(ingredients,columns=['ID','Peptide'])
        print("Homebrow is Done!")
    recipie[['db', 'UniprotID','ID2']] = recipie['ID'].str.split('|', 2, expand=True)
    recipie[['Gene','Identification']] = recipie['ID2'].str.split('_', 1, expand=True)
    recipie.drop(columns=['ID', 'ID2',"db"], inplace=True)
    return(recipie)