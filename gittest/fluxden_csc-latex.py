#%%
import pandas as pd
pd.set_option('display.max_columns', None) #显示所有行
pd.set_option('display.max_rows', None)
df_o = pd.read_csv("../fluxdensity/IntFlux_final.csv", header=None)

df_n = df_o.transpose()
# df_filter = df_tran[df_tran.loc[3] == "spw.mask.pwn.fits"]
df_n.iloc[4:7] *= 1e3

df_latex = df_n.to_latex(index=False)
df_n.head(100)
# print(df_latex)

# make a new df
# df_latex = pd.DataFrame(index=["central frequency (MHz)",
#                                "Whole flux density (mJy)",
#                                "Whole flux density error (mJy) ",
#                                "Surface brightness at 1 GHz",
#                                "Shell flux density (mJy)",
#                                "Shell flux density error (mJy)",
#                                "PWN flux density (mJy)",
#                                "PWN flux density error (mJy)"
#                                ],
#                         columns=["MFS", "spw1", "spw2", "spw3", "spw4"])


# df_format = df_latex.copy()
# df_format.loc[-1] = ['\hline'] * len(df_latex)
# df_format.sort_index(inplace=True)

# df_latex = df_format.to_latex(column_format="lccc")
# print(repr(df_latex))
# %%
