#%%
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from globalSetting import *


set_plot_settings()


# 读取数据，分隔符为逗号，跳过第二行，
df = pd.read_csv("../Data/RM-DM_in10deg.txt", sep=',', skiprows=[1])

# filtered RAD<2
df = df[df["RAD"] < 2]

# 选取J1400的数据
df_J1400 = df[df["NAME"] == "J1400-6325"]

RM_J1400 = -670
DM_J1400 = df_J1400["DM"].values[0]
dist_J1400 = df_J1400["DIST"].values[0]

df_RM = df.dropna(subset=['RM', 'DIST'])
df_DM = df.dropna(subset=['DM', 'DIST'])

# df_J1400_RM = df_RM[df_RM["NAME"] == "J1400-6325"]
# df_J1400_DM = df_DM[df_DM["NAME"] == "J1400-6325"]

# 提取RM、Dist和DM列的数据
RM = df_RM["RM"]
dist_RM = df_RM['DIST']


fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(211)

ax1.scatter(dist_RM, RM, color='royalblue', marker='o', s=100, alpha=0.5)

# 绘制J1400的RM和DM的散点图，标记为红色
# ax1.scatter(dist_J1400, RM_J1400, color='red', marker='+', s=100, linewidth = 3, label='PSR J1400-6325')
ax1.axhline(y=RM_J1400, color='tomato', linewidth = 3, linestyle='--', label='SNR G310.6-1.6')

ax1.set_ylabel(r"RM ($\rm{rad~m ^{-2}})$")
# ax1.set_title('RM vs Distance')
ax1.legend()


# 提取DM列的数据
DM = df_DM["DM"]
dist_DM = df_DM['DIST']
ax2 = fig.add_subplot(212)

# 绘制DM vs Distance的散点图
ax2.scatter(dist_DM, DM, color='royalblue', marker='o', s=100, alpha=0.5)
ax2.scatter(dist_J1400, DM_J1400, color='red', marker='+', s=100, linewidth = 3,label='PSR J1400-6325',)

ax2.set_xlabel('Distance (kpc)' )
ax2.set_ylabel(r"DM ($\rm{cm^{-3}~pc})$")
# ax2.set_title('DM vs Distance')

ax2.legend()

df_RMDM= df.dropna(subset=['RM', 'DM'])
RM = df_RMDM["RM"]
DM = df_RMDM["DM"]

# ax3 = fig.add_subplot(213)


# ax3.scatter(RM, DM, color='royalblue', marker='o', s=100, alpha=0.5)
# ax3.scatter(RM_J1400, DM_J1400, color='red', marker='+', s=100, linewidth = 3, alpha=0.8, label='PSR J1400-6325' )
# ax3.set_xlabel(r"RM ($\rm{rad~m ^{-2}})$")
# ax3.set_ylabel(r"DM ($\rm{cm^{-3}~pc})$")
# ax3.legend()


plt.subplots_adjust(wspace =0, hspace = 0.3)
# sv_fig("../figures/RM-DM-dist")
plt.show()



# %%
