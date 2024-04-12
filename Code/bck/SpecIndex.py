# %%

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import pandas as pd
from matplotlib.ticker import ScalarFormatter
import inspect
def sv_fig(name):
    sv = plt.savefig(name + '.png', bbox_inches='tight', dpi = 300)
    sv = plt.savefig(name + '.pdf', bbox_inches='tight')
    return sv


config = {
    "font.family":'serif',
    "font.size": 20,
    "mathtext.fontset":'stix',
    "font.serif": ['SimSun'], # simsun字体中文版就是宋体
}
plt.rcParams['font.size'] = 20
plt.rc('font',family='Times New Roman')


def show_specindex(data, plot=True):
    freq = data[:, 0]
    flux = data[:, 1]
    err  = data[:, 2]


    freq_lg = np.log10(freq)
    flux_lg = np.log10(flux)
    err_u = np.log10(flux + err) - np.log10(flux)
    err_l = abs(np.log10(flux - err) - np.log10(flux))
    err_lg = np.maximum(err_u, err_l) # 取两个数组中的最大值
    print(err_u,"\n", err_l, "\n", err_lg)

    model = models.Linear1D()
    fitter = fitting.LinearLSQFitter()
    best_fit = fitter(model, freq_lg, flux_lg, weights=1.0/err_lg) # weight 加权重，将较小的权重加给误差更大的数据点 **2
    
    specindex = best_fit.parameters[0]
    b = best_fit.parameters[1]
    print("specindex",specindex,"b", b)


    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(freq_lg, flux_lg, yerr=err_lg, fmt='o', ecolor='skyblue', color='darkblue', elinewidth=4, capsize=10, ms=8)
        

        # Plotting the best fit line
        y_fit = best_fit(freq_lg)#10 ** best_fit(np.log10(freq))
        ax.plot(freq_lg, y_fit, color='orange', linewidth=3, linestyle = "--", label=f'Spectral index: {specindex:.2f}')

        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Flux (mJy)')
        # ax.set_title('Spectral Index')

        # 设置横坐标和纵坐标为对数刻度
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.legend(loc='best')

    return ax


data_PWN = np.array([
    # (0.1,3.27915042,0.003),
    (8.354910000000E+08, 0.2540, 0.0039),
    (9.074910000000E+08, 0.2467, 0.0030),
    (9.794910000000E+08, 0.2392, 0.0024),
    (1.051490000000E+09, 0.2296, 0.0023)
])
ax = show_specindex(data_PWN, True)



# data_SNR = np.array([
#     (8.354910000000E+08, 0.0365, 0.0052),
#     (9.074910000000E+08, 0.0355, 0.0039),
#     (9.794910000000E+08, 0.0334, 0.0032),
#     (1.051490000000E+09, 0.0311, 0.0030)
# ])

# ax = show_specindex(data_SNR, True)




data_shell  = np.array([
    (8.354910000000E+08, 0.0425, 0.0053),
    (9.074910000000E+08, 0.0413, 0.0042),
    (9.794910000000E+08, 0.0390, 0.0036),
    (1.051490000000E+09, 0.0362, 0.0032)
])

ax = show_specindex(data_shell, True)

# ax.errorbar(data_orig[2,0], data_orig[2,1], yerr=data_orig[2,2], fmt='o', ecolor='pink', color='darkred', elinewidth=4, capsize=10, ms=8, label = "MFS")

# sv_fig("../Output/SpecIndexFig_SNR" )
plt.show()



# %%
