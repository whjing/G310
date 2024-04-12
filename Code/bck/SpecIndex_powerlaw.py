# %%
##################
## Note
# 还没弄好，不知道这个要怎么改。
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import pandas as pd
from matplotlib.ticker import ScalarFormatter
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


#PWN perfect
data_orig = np.array([
    # (0.1,3.27915042,0.003),
    (8.354910000000E+08, 0.2540, 0.0039),
    (9.074910000000E+08, 0.2467, 0.0030),
    (9.794910000000E+08, 0.2392, 0.0024),
    (1.051490000000E+09, 0.2296, 0.0023)
])





def show_specindex(data, plot=True):
    freq = data[:, 0]
    flux = data[:, 1]
    err = data[:, 2]


    freq_lg = np.log10(freq)
    flux_lg = np.log10(flux)
    err_u = np.log10(flux + err) - np.log10(flux)
    err_l = abs(np.log10(flux - err) - np.log10(flux))
    err_lg = np.maximum(err_u, err_l) # 取两个数组中的最大值
    print(err_u,"\n", err_l, "\n", err_lg)

    model = models.PowerLaw1D()
    fitter = fitting.SimplexLSQFitter()
    best_fit = fitter(model, freq, flux, weights=1.0/err**2) # weight 加权重，将较小的权重加给误差更大的数据点
    print(best_fit)
    specindex = best_fit.parameters[0]


    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(freq, flux, yerr=err_lg, fmt='o', ecolor='skyblue', color='darkblue', elinewidth=4, capsize=10, ms=8)
        

        # Plotting the best fit line
        y_fit = best_fit(flux)
        y_fit = 10 ** best_fit(np.log10(freq))
        ax.plot(freq, y_fit, color='orange', linewidth=3, linestyle = "--", label=f'Spectral index: {specindex:.2f}')

        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Flux (mJy)')
        # ax.set_title('Spectral Index')

        # 设置横坐标和纵坐标为对数刻度
        ax.set_xscale('log')
        ax.set_yscale('log')

        
        
        # plt.grid()

        

    return ax

# data = data_orig[[0,3]]
# data_point = (9.43E+08, 2.58, 0.002)
data = data_orig
ax = show_specindex(data, True)
# ax.errorbar(data_orig[2,0], data_orig[2,1], yerr=data_orig[2,2], fmt='o', ecolor='pink', color='darkred', elinewidth=4, capsize=10, ms=8, label = "MFS")

ax.legend(loc='best') 
sv_fig("../Output/SpecIndexFig_SNR" )
plt.show()



# %%
