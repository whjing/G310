#%%
import numpy as np
import pandas as pd
#%%

greencat = pd.read_csv("~/academic/Cat/greencat.csv")
# set "name" as index
greencat.set_index(["name"], inplace=False)


xssnr = pd.read_csv("../cat/XSSNR.csv")

xssnr_n = pd.merge(xssnr,greencat, left_on="Gname", right_on="name", how = "inner")
xssnr_n = xssnr_n[["name_x", "Gname", "RA", "Dec", "size", "type", 'flux',"s"]]
# xssnr_n.rename(columns={"flux": "fluxDensity_radio_1GHz"}, inplace=True)
# xssnr_n = pd.concat([xssnr,greencat], axis = 1)#,join="inner")
xssnr_n.head(n = 8)

# xssnr_n.to_csv("../cat/XSSNR.csv",mode='w', index="name_x")

# %%
import numpy as np
import pandas as pd


xssnr = pd.read_csv("../cat/XSSNR_old.csv", skiprows=[1])

file_snrcat = "/Users/jing/academic/Cat/SNRcat.csv"
snrcat = pd.read_csv(file_snrcat, sep=";", skiprows=2)
snrcat.info()
snrcat["name"] = snrcat.apply(lambda row: f'G{float(row["G"][1:6])}{row["G"][6]}{float(row["G"][7:11])}', axis=1)

xssnr_n = pd.merge(xssnr,snrcat, left_on="Gname", right_on="name", how = "inner")

# snrcat.loc[:, "name"] = snrcat["G"].str[1:].replace('.0', '.')

# text = open(file_snrcat).readlines()[0:5]
# print(text)


xssnr_n.drop(columns=["G", "J2000_ra (hh:mm:ss)", "J2000_dec (dd:mm:ss)", "J2000_from", "id_uncertain", "id_new","size_imprecise"], inplace=True)

xssnr_n[["name_x","size", "size_radio","size_X"]].head(10)

# xssnr_n.to_csv("../cat/XSSNR_n.csv",mode='w', index="name_x")
# %%
def mergerCat():
    file_snrcat = "/Users/jing/academic/Cat/SNRcat.csv"
    file_greencat = "/Users/jing/academic/Cat/greencat.csv"
    snrcat = pd.read_csv(file_snrcat, sep=";", skiprows=2)
    greencat = pd.read_csv(file_greencat)
    snrcat["name"] = snrcat.apply(lambda row: f'G{float(row["G"][1:6])}{row["G"][6]}{float(row["G"][7:11])}', axis=1)
    catMerged = pd.merge(snrcat, greencat, left_on="name", right_on="name", how = "inner", suffixes=("_HE", "_G"))
    print("ss", greencat["type"].unique())
    return catMerged


allsnr = mergerCat()
all_pler = allsnr["type_HE"].str.contains("plerionic")
all_C = allsnr["type_G"].str.contains("C") 
all_pwn = allsnr["context"].str.contains("PWN")
all_snr = allsnr["context"].str.contains("SNR")

# HEcat or Gcat
# CCsnr_large = allsnr[all_pler  all_C]

# HEcat and Gcat 且 比较明显的SNR或PWN
CCsnr = allsnr[all_pler | all_C ]#& (all_pwn|all_snr)]
# CCsnr.loc[CCsnr["size_radio"].notna(), ]

# cat = CCsnr.loc[CCsnr["size_radio"].notna(), ]
# cat["size_radio"]
# cat
pd.options.display.max_columns = None

# 所以其实size_radio 并不是实际radio所探测到的大小，而是从Green cat中直接抄来的
print(CCsnr.columns)
# snrcat["name"] = snrcat.apply(lambda row: f'G{float(row["G"][1:6])}{row["G"][6]}{float(row["G"][7:11])}', axis=1)

CCsnr["size_left"] = CCsnr.apply(lambda row: f"{float(row['size'].replace('?','').split('x')[0])}", axis=1)
CCsnr["size_right"] = CCsnr.apply(lambda row: f"{float(row['size'].replace('?','').split('x')[1]) if len(row['size'].replace('?','').split('x')) > 1 else float(row['size'].replace('?','').split('x')[0])}", axis=1)

# CCsnr["size_float"] = CCsnr["size"].str.split(".")[0]
CCsnr = CCsnr [["name", "RA", "Dec", "size", "size_left","size_right","size_radio", "size_X", "type_HE","type_G", 'flux',"s",'age_min (yr)', 'age_max (yr)', 'distance_min (kpc)','context', ]]
cc1 = CCsnr.sort_values(by="size_left", ascending=True)
# cc1[cc1["name"].str.contains("G21")]

cc1
# print(cc1.to_latex())

# print(CCsnr_large.shape)
# print(CCsnr.shape)

# CCsnr.head(10)

# 所以筛选了一大堆之后，目前混合型超新星遗迹中，只有G310由于其小尺寸，所以没有探测到

# %%
