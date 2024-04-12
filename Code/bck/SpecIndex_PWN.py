#%%
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

data_orig = np.array([
    (0.843, 217.4, 9.4),
    (0.944, 242, 5.6),
    (2.4  , 600  , 0 ),
    (4.5  , 113  , 10),
    (4.5  , 131  , 10),
    (5.5  , 150  ,  9),
    (9.   , 147  ,  9)
])

data_list = [
    "SUMSS",
    "EMU",
    "Parkes",
    "PMN",
    "PMN",
    "ATCA",
    "ATCA",
]

data = data_orig
df = pd.DataFrame(data, columns=['Freq (GHz)', 'Flux (mJy)', 'F_err (mJy)'])
df['Obs'] = data_list
print(df)
#%%
def show_specindex(data, plot=True):
    freq_orig = data[:, 0]
    flux_orig = data[:, 1]
    error = data[:, 2]

    freq = np.log10(freq_orig)
    flux = np.log10(flux_orig)

    model = models.Linear1D()
    fitter = fitting.LinearLSQFitter()
    best_fit = fitter(model, freq, flux)
    specindex = best_fit.parameters[0]

    df = pd.DataFrame(data, columns=['Freq (GHz)', 'Flux (mJy)', 'Error'])

    # 打印DataFrame
    print(df)

    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(freq_orig, flux_orig, yerr=error, fmt='o', ecolor='skyblue', color='darkblue', elinewidth=4, capsize=10, ms=8)
        

        # Plotting the best fit line
        x_fit = np.linspace(min(freq_orig), max(freq_orig), 100)
        y_fit = 10 ** best_fit(np.log10(x_fit))
        ax.plot(x_fit, y_fit, color='orange', linewidth=3, linestyle = "--", label=f'Spectral index: {specindex:.2f}')

        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Flux (mJy)')
        # ax.set_title('Spectral Index')

        # 设置横坐标和纵坐标为对数刻度
        ax.set_xscale('log')
        ax.set_yscale('log')

        
        
        # plt.grid()

        

    return ax

data = data_orig[[0,3]]

# data = data_orig[[0,4,5,6]]
ax = show_specindex(data, True)
ax.errorbar(data_orig[1,0], data_orig[1,1], yerr=data_orig[1,2], fmt='o', ecolor='pink', color='darkred', elinewidth=4, capsize=10, ms=8, label = "EMU")

ax.legend(loc='best') 
# sv_fig("../Output/SpecIndex_PWN_all" )
# data = data[[0,1,4,5,6,]]
# show_specindex(data, True)
# sv_fig("../Output/SpecIndex_PWN_all_my" )
# sv_fig("../Output/SpecIndex_PWN_all" )
plt.show()

# %%


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

# data_orig = np.array([
#  (0.872, 261.2 , 8),
#  (1.016, 183.2 , 7)
# ])

# PWN
# data_orig = np.array([
#     (8.354910000000E+08, 0.2280, 0.0033),
#     (9.074910000000E+08, 0.2132, 0.0028),
#     (9.794910000000E+08, 0.2005, 0.0023),
#     (1.051490000000E+09, 0.1879, 0.0022)
# ]

# #PWN
# data_orig = np.array([
#     (8.354910000000E+08,  0.2160 ,  0.0034),
#     (9.074910000000E+08,  0.2046 ,  0.0025),
#     (9.794910000000E+08,  0.1920 ,  0.0019),
#     (1.051490000000E+09,  0.1771 ,  0.0018)
# ])


# #PWN
# data_orig = np.array([
#     (8.354910000000E+08,  1.850561940236e-1  ,  0.0034),
#     (9.074910000000E+08,  1.757892517605e-1 ,  0.0025),
#     (9.794910000000E+08,  1.644564998134e-1 ,  0.0019),
#     (1.051490000000E+09,  1.330880247135e-1  ,  0.0018)
# ])

# #PWN calIntFlux t
# data_orig = np.array([
#     (8.354910000000E+08,  0.2584, 0.0088),
#     (9.074910000000E+08,  0.2345, 0.0075),
#     (9.794910000000E+08,  0.2132, 0.0071),
#     (1.051490000000E+09,  0.1920, 0.0063)
# ])

# #PWN carta core
# data_orig = np.array([
#     (8.354910000000E+08,  2.035118459413e-1 , 0.0088),
#     (9.074910000000E+08,  1.849918721423e-1 , 0.0075),
#     (9.794910000000E+08,  1.673754875633e-1, 0.0071),
#     (1.051490000000E+09,  1.498231476675e-1 , 0.0063)
# ])



#PWN calIntFlux t
# data_orig = np.array([
#     (8.354910000000E+08,  1.552247005514e-1, 0.0088),
#     (9.074910000000E+08,  1.415733314748e-1, 0.0075),
#     (9.434910000000E+08,  1.354742741898e-1, 0.0020),
#     (9.794910000000E+08,  1.296215663308e-1, 0.0071),
#     (1.051490000000E+09,  1.168171771980e-1, 0.0063)
# ])

# data_orig = np.array([
#     (8.354910000000E+08,  1.552247005514e-1, 0.0088),
#     (9.074910000000E+08,  1.415733314748e-1, 0.0075),
#     (9.434910000000E+08,  1.354742741898e-1, 0.0020),
#     (9.794910000000E+08,  1.296215663308e-1, 0.0071),
#     (1.051490000000E+09,  1.168171771980e-1, 0.0063)
# ])


# data_orig = np.array([
#     (8.354910000000E+08,  	1.109695771162e-1 , 0.0088),
#     (9.074910000000E+08,  	1.017408149495e-1 , 0.0075),
#     (9.794910000000E+08,  	9.336624863495e-2 , 0.0071),
#     (1.051490000000E+09,  	8.473634132353e-2, 0.0063)
# ])


# data_orig = np.array([
#     (8.354910000000E+08,  	1.331251161229e-1 , 0.0088),
#     (9.074910000000E+08,  	1.300889489423e-1 , 0.0075),
#     (9.794910000000E+08,  	1.270381333028e-1  , 0.0071),
#     (1.051490000000E+09,  	1.228648768901e-1 , 0.0063)
# ])


# data_orig = np.array([
#     (8.354910000000E+08,  	1.360247171492e-2 , 0.0088),
#     (9.074910000000E+08,  	1.401424221098e-2 , 0.0075),
#     (9.794910000000E+08,  	1.272597251205e-2  , 0.0071),
#     (1.051490000000E+09,  	1.055372405233e-2 , 0.0063)
# ])



# data_orig = np.array([
#     (8.354910000000E+08,  		1.236076080322e-2, 0.0088),
#     (9.074910000000E+08,  	1.103884079130e-2  , 0.0075),
#     (9.794910000000E+08,  		9.988116587718e-3 , 0.0071),
#     (1.051490000000E+09,  8.455366475255e-3 , 0.0063)
# ])

#PWN perfect
# data_orig = np.array([
#     # (0.1,3.27915042,0.003),
#     (8.354910000000E+08, 0.2540, 0.0039),
#     (9.074910000000E+08, 0.2467, 0.0030),
#     (9.794910000000E+08, 0.2392, 0.0024),
#     (1.051490000000E+09, 0.2296, 0.0023)
# ])


data_orig  = np.array([
    (8.354910000000E+08, 0.035420, 0.0053),
    (9.074910000000E+08, 0.034362, 0.0042),
    (9.794910000000E+08, 0.032430, 0.0036),
    (1.051490000000E+09, 0.030176, 0.0032)
])

# # ERS 1

# data_orig = np.array([
#     (8.354910000000E+08, 0.8328, 0.0016),
#     (9.074910000000E+08, 0.7878, 0.0012),
#     (9.794910000000E+08, 0.7524, 0.0011),
#     (1.051490000000E+09, 0.7147, 0.0009)
# ])


# # ERS 2
# data_orig = np.array([
#     (8.354910000000E+08, 0.0094, 0.0011),
#     (9.074910000000E+08, 0.0088, 0.0009),
#     (9.794910000000E+08, 0.0087, 0.0008),
#     (1.051490000000E+09, 0.0084, 0.0007)
# ])

# # ERS 3

# data_orig = np.array([
#     (8.354910000000E+08, 0.0026, 0.0009),
#     (9.074910000000E+08, 0.0027, 0.0007),
#     (9.794910000000E+08, 0.0026, 0.0006),
#     (1.051490000000E+09, 0.0026, 0.0005)
# ])






# SNR great

# data_orig = np.array([
#     (8.354910000000E+08, 0.0365, 0.0052),
#     (9.074910000000E+08, 0.0355, 0.0039),
#     (9.794910000000E+08, 0.0334, 0.0032),
#     (1.051490000000E+09, 0.0311, 0.0030)
# ])






# SNR
# data_orig = np.array([
#     # (8.354910000000E+08, 0.0284, 0.0040),
#     (9.074910000000E+08, 0.0264, 0.0033),
#     (9.794910000000E+08, 0.0249, 0.0028),
#     (1.051490000000E+09, 0.0232, 0.0026)
# ])





# # PSR
# data_orig = np.array([
#     (1.374, 0.25, 0),
#     (3.1,   0.11, 0),
#     (4,     0.08, 0),
# ])


def show_specindex(data, plot=True):
    freq_orig = data[:, 0]
    flux_orig = data[:, 1]
    error = data[:, 2]

    freq = np.log10(freq_orig)
    flux = np.log10(flux_orig)

    model = models.Linear1D()
    fitter = fitting.LinearLSQFitter()
    best_fit = fitter(model, freq, flux)
    specindex = best_fit.parameters[0]
    b = best_fit.parameters[1]
    print("specindex",specindex,"b", b)
    freq_1_orig = 1e9
    flux_1_orig = 10 ** best_fit(np.log10(freq_1_orig))
    print(flux_1_orig)



    # df = pd.DataFrame(data, columns=['Freq (GHz)', 'Flux (mJy)', 'Error'])

    # 打印DataFrame
    # print(df)

    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.errorbar(freq_orig, flux_orig, yerr=error, fmt='o', ecolor='skyblue', color='darkblue', elinewidth=4, capsize=10, ms=8)
        

        # Plotting the best fit line
        x_fit = np.linspace(min(freq_orig), max(freq_orig), 100)
        y_fit = 10 ** best_fit(np.log10(x_fit))
        ax.plot(x_fit, y_fit, color='orange', linewidth=3, linestyle = "--", label=f'Spectral index: {specindex:.2f}')

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
# sv_fig("../Output/SpecIndex_PWN_all" )
# data = data[[0,1,4,5,6,]]
# show_specindex(data, True)
# sv_fig("../Output/SpecIndex_PWN_all_my" )
# sv_fig("../Output/SpecIndexFig_SNR" )

plt.show()



# %%
print((0.0345 + 0.0425) * 0.46)
print((0.0334 + 0.0413) * 0.46)
print((0.0315 + 0.0390) * 0.46)
print((0.0294 + 0.0362) * 0.46)
# %%
