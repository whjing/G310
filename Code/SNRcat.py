#%%
import pandas as pd
file_snrcat = "~/academic/Cat/SNRcat.csv"
import numpy as np

snrcat = pd.read_csv(file_snrcat, delimiter = ";", skiprows=2)
pd.set_option('display.max_columns', None)
# print(snrcat)


snrcat_plerionic = snrcat["type"].str.contains("thermal") 
snrcat_sizeX = snrcat["size_X"].notna()
snrcat_sizeradio = snrcat["size_radio"].notna()

snrcat_select = snrcat[snrcat_plerionic ]#& snrcat_sizeX & snrcat_sizeradio]
# snrcat_select.to_csv("../data/snrcat_plerionic.csv",sep=",")
snrcat_select.info()
# print([snrcat[snrcat["G"]=="G310.6-01.6"]])
snrcat_select.head()
# %%
