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

def MeatWrapper(workbook,old,new,output):
    meat = openpyxl.load_workbook(workbook)
    wrap = meat[old]
    wrap.title = new
    meat.save(output)

def Excel_Mapper(list_dfs, xls_path):
    with ExcelWriter(xls_path) as writer:
        for n, df in enumerate(list_dfs):
            df.to_excel(writer,'sheet%s' % n)
        writer.save()

def PEAKS_Importer(csv,drop_OG=True):
    df=pd.read_csv(csv)
    df[["Protein","Y"]]=df['Protein Accession'].str.split("|",1,expand=True)
    df[["Gene","Species"]]=df['Y'].str.split("_",1,expand=True)
    if drop_OG==True:
        df.drop(columns=['Y', 'Protein Accession',"Found By"], inplace=True)
    else:
        df.drop(columns=['Y',"Found By"], inplace=True)
    return(df)

#Handles Up to 3 Replicates (t_id1-3) per Selection, used for removing amino acids from N-terminal during PEAKS exports with enzyme that cleave at C-terminal. 
# Can be used to remove M from N-terminal of peptides produced by enzymes which cleave at N-terminal of target aa. 
def Butcher(df,ident1=None,ident2=None,ident3=None,t_id1=None,t_id2=None,t_id3=None,t_value=0,acid=["J","Z"],labels=list(),excel_mapper=True,excel_name=None):
    raw=df.loc[:,df.columns.str.contains(ident1)]
    tag=df[labels]
    raw = pd.concat([raw, tag], axis=1)
    raw["Peptide"]= raw["Peptide"].str.replace('\W+',"")
    raw["Peptide"]= raw["Peptide"].str.replace('\d+',"")
    raw["Peptide"]= raw["Peptide"].apply(lambda x : x[1:] if x.startswith(tuple(acid)) else x)
    cuts=raw.loc[:,raw.columns.str.contains(ident2)]
    if ident3 != None:
        cuts=cuts.loc[:,cuts.columns.str.contains(ident3)]
    cuts = pd.concat([cuts, tag], axis=1)
    blade = cuts.filter(regex=r'^AREA').isin(['0']).all(axis=1)
    cuts=cuts.loc[~blade]
    cuts.reset_index(inplace=True)
    excels=[raw,cuts]
    Excel_Mapper(excels,excel_name +".xlsx")
    return raw, cuts


def Peptide_Mass(peptide):
    mass = {  "A": 71.037114,
    "R": 156.101111,
    "N": 114.042927,
    "D": 115.026943,
    "C": 103.009185,
    "Q": 129.042593,
    "E": 128.058578,
    "G": 57.021464,
    "H": 137.058912,
    "I": 113.084064,
    "L": 113.084064,
    "K": 128.094963,
    "M": 131.040485,
    "F": 147.068414,
    "P": 97.052764,
    "S": 87.032028,
    "T": 101.047679,
    "W": 186.079313,
    "Y": 163.06332,
    "V": 99.068414,
    }
    mass_list = [mass.get(aa,0.0)for aa in peptide]
    pep_mass=sum(mass_list)
    return pep_mass 

def ButcherShop(df, target,identifier, rule, min_length=7,exception=None,max_length=100, pH=2.0, min_charge=2,missed=0):
    raw=df[[target, identifier]].set_index(identifier).to_dict()[target]
    string_catcher=re.compile(r'^([A-Z]+)$')
    pep_dict = {}
    pep_dict_list = []
    print(f'You order is being processed and the butcher is preparing your {rule}-cut protein(s)!')
    print("The butcher is working...")
    for gene,peptide in raw.items():
        pep_dict[gene] = list(parser.cleave(peptide,rule=rule,min_length=min_length,exception=exception,missed_cleavages=missed))
    for k, lst in pep_dict.items():    
        d = {}
        for i in range(len(lst)):
            d.update({k: lst[i]})
            d.update({'gene':k})
            d.update({'aa_comp': dict(parser.amino_acid_composition(lst[i]))})
            d.update({'peptide': re.findall(string_catcher,lst[i])})
            d.update({'Length': len(lst[i])})
            d.update({'z': int(round(electrochem.charge(lst[i], pH=pH)))})
            d.update({'Mass': int(Peptide_Mass(lst[i]))})
            if d["z"] > 0:
                d.update({'m/z': d["Mass"]/d['z']})
            new_d = d.copy()
            pep_dict_list.append(new_d)
    print("Trimming the cuts....")
    pep_dict_list = [peptide for peptide in pep_dict_list if peptide['Length'] <= int(max_length)]
    print("Weighing the cuts...")
    pep_dict_list = [peptide for peptide in pep_dict_list if peptide['z'] >= int(min_charge)]
    print(f'Order is up! You have acquired {len(pep_dict_list)} peptides that are between {min_length} and {max_length} amino acids!')
    return pep_dict_list


def DeliShop(z,meat_package=False):
    #z = list of dictionaries, keys must be equal thus will drop keys which are not cosistent between dictionaries
    # use after ButcherShop
    #returns dataframe
    key_intersect = set(z[0].keys()).intersection(set(z[1].keys()))
    zz = [{key:value for (key,value) in dicts.items() if key in key_intersect} for dicts in z]
    ham = pd.DataFrame(zz)
    if meat_package == True:
        ham_counts=ham.groupby('gene').size().reset_index(name='counts')
        df=df.merge(ham_counts,how='left', on=['gene'])
    return df








