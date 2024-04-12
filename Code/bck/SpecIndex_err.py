
# %%


import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import curve_fit



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

plt.rc('font', family='Times New Roman')

def linear_func(x, a, b):
    return a * x + b


def show_specindex(data, plot=True):
    freq_orig = data[:, 0]
    flux_orig = data[:, 1]
    flux_error = data[:, 2]

    # print(f"flux_up {flux_up} and flux_low {flux_low}")

    flux_lowmax  = np.log10(data[0, 1]  + data[0, 2] )
    flux_lowmin  = np.log10(data[0, 1]  - data[0, 2] )
    flux_highmax = np.log10(data[-1, 1] + data[-1, 2])
    flux_highmin = np.log10(data[-1, 1] - data[-1, 2])
    print(
        flux_lowmax ,
        flux_lowmin ,
        flux_highmax,
        flux_highmin,

    )



    # print(f"flux max min: {flux_maxmin}")



    freq = np.log10(freq_orig)
    flux = np.log10(flux_orig)
    print(flux)

    flux_up  = flux_orig + flux_error
    flux_low = flux_orig - flux_error
    specIndex_max = np.diff(
        np.array([
        (freq[0], flux_lowmax),
        (freq[-1], flux_highmin),
        ])
        )
    specIndex_min = flux_lowmin - flux_highmax / ( freq[0] - freq[-1])
    print(f"specIndex max {specIndex_max}")
    print(f"specIndex min {specIndex_min}")
    # flux_error = abs(np.log10(flux_error))

    # 使用 curve_fit 进行拟合
    popt, pcov = curve_fit(linear_func, freq, flux , sigma=flux_error)
    specindex, b_fit = popt
    specindex_err = np.float64(np.sqrt(np.diag(pcov))[0])


    # # 创建一个 Linear1D 模型
    # model = models.Linear1D()

    # # 创建一个 LinearLSQFitter，并传入权重参数
    # fitter = fitting.LinearLSQFitter()

    # # 使用 flux_err 作为权重进行拟合
    # best_fit = fitter(model, freq, flux, weights=1/flux_error)

    # # 获取拟合结果中的斜率参数
    # specindex = best_fit.parameters[0]
    # specindex_err = np.sqrt(np.diag(fitter.fit_info["residuals"]))[0]


    # model = models.Linear1D()
    # fitter = fitting.LinearLSQFitter()
    # best_fit = fitter(model, freq, flux)
    # specindex = best_fit.parameters[0]

    print(f"Spectral index: {specindex:.2f} flux_error {specindex_err:.2f}")

    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(freq_orig, flux_orig, yerr=flux_error, fmt='o', ecolor='skyblue', color='darkblue', elinewidth=4, capsize=10, ms=8)
        ax.fill_between(freq_orig, flux_low, flux_up, facecolor='skyblue', alpha = 0.5)
        
        # Plotting the best fit line
        x_fit = np.linspace(min(freq_orig), max(freq_orig), 100)
        # y_fit = 10 ** best_fit(np.log10(x_fit))
        y_fit = 10 ** (specindex * np.log10(x_fit) + b_fit)
        ax.plot(x_fit, y_fit, color='orange', linewidth=3, linestyle = "--", label = f"Spectral index: ${specindex:.2f} \pm {specindex_err:.2f}$")
        

        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Flux (mJy)')
        # ax.set_title('Spectral Index')

        # 设置横坐标和纵坐标为对数刻度
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        
        # plt.grid()

        

    return ax


# #PWN perfect
# data_orig = np.array([
#     (8.354910000000E+08, 0.2540, 0.0039),
#     (9.074910000000E+08, 0.2467, 0.0030),
#     (9.794910000000E+08, 0.2392, 0.0024),
#     (1.051490000000E+09, 0.2296, 0.0023)
# ])

# SNR great

# data_orig = np.array([
#     (8.354910000000E+08, 0.0365, 0.0052),
#     (9.074910000000E+08, 0.0355, 0.0039),
#     (9.794910000000E+08, 0.0334, 0.0032),
#     (1.051490000000E+09, 0.0311, 0.0030)
# ])

# SNR Carta
# data_orig = np.array([
#     (8.354910000000E+08, 0.05446184990980002 , 1e-5),
#     (9.074910000000E+08, 0.0538466529919     , 1e-5),
#     (9.794910000000E+08, 0.04759179779739997 , 1e-5),
#     (1.051490000000E+09, 0.039860834681099994, 1e-5),
# ])


data_orig = np.array([
    (8.354910000000E+08,2.6962285593600015, 1e-5),
    (9.074910000000E+08,2.661314559340001, 1e-5),
    (9.794910000000E+08,2.3654786462399997, 1e-5),
    (1.051490000000E+09,1.998672611089999, 1e-5),
])

# SNR with bkg
data_orig = np.array([
    (8.354910000000E+08,0.0387,0.0050),
    (9.074910000000E+08,0.0380,0.0039),
    (9.794910000000E+08,0.0334,0.0034),
    (1.051490000000E+09,0.0275,0.0030),
])



# 8.35491E+08
# 9.07491E+08
# 9.79491E+08
# 1.05149E+09

# data = data_orig[[0,3]]
# data_point = (9.43E+08, 2.58, 0.002)
data = data_orig
ax = show_specindex(data, True)
# ax.errorbar(data_orig[2,0], data_orig[2,1], yerr=data_orig[2,2], fmt='o', ecolor='pink', color='darkred', elinewidth=4, capsize=10, ms=8, label = "MFS")

ax.legend(loc='best') 
# sv_fig("../Output/SpecIndex_PWN_all" )
# data = data[[0,1,4,5,6,]]
# show_specindex(data, True)
# sv_fig("../Output/SpecIndex_PWN_all_my" )
# sv_fig("../Output/SpecIndexFig_SNR" )
plt.show()



# %%
