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


def show_specindex(data, plot=True):
    data = data[~np.any(data == 0, axis=1)]
    freq_orig = data[:, 0]
    flux_orig = data[:, 1]
    
    print(data)

    freq = np.log10(freq_orig)
    flux = np.log10(flux_orig)

    model = models.Linear1D()
    fitter = fitting.LinearLSQFitter()
    best_fit = fitter(model, freq, flux)
    specindex = best_fit.parameters[0]


    if plot:
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.scatter(freq_orig, flux_orig,)

        if len(data[0,:]) == 3:
            error = data[:, 2]
            ax.errorbar(freq_orig, flux_orig, yerr=error, fmt='o', ecolor='skyblue', color='darkblue', elinewidth=4, capsize=10, ms=8)

        # Plotting the best fit line
        x_fit = np.linspace(min(freq_orig), max(freq_orig), 100)
        y_fit = 10 ** best_fit(np.log10(x_fit))
        ax.plot(x_fit, y_fit, color='orange', linewidth=3, linestyle = "--", label=f'Spectral index: {specindex:.2f}')

        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Flux (mJy)')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc='best') 

    return ax

# data_convolved = np.loadtxt("../Data/icube.convolved.single-pix-spectra.txt")

# show_specindex(data_convolved, True)

# # data_unconvolved = np.loadtxt("../Data/icube.unconvolved.single-pix-spectra.txt")
# # show_specindex(data_unconvolved, True)

# sv_fig("../Output/1")
# data_convolved_flag = data_convolved[40:,:]

# show_specindex(data_convolved_flag, True)

# # data_unconvolved_flag = data_unconvolved[40:,:]
# # show_specindex(data_unconvolved_flag, True)
# sv_fig("../Output/2")


# file_PWN = "../Data/i.cube.convolved.PWN-spectra.txt"
# data_PWN = np.loadtxt(file_PWN)
# show_specindex(data_PWN, True)
# sv_fig("../Output/SpecIndex_PWN_orig")
# data_PWN_flag = data_PWN[55:]
# show_specindex(data_PWN_flag, True)
# sv_fig("../Output/SpecIndex_PWN_flag")

# file_PWNcore = "../DAta/i.cube.convolved.PWN-core-spectra.txt"
# data_PWN = np.loadtxt(file_PWNcore)
# show_specindex(data_PWN, True)
# sv_fig("../Output/SpecIndex_PWN-core_orig")
# data_PWN_flag = data_PWN[55:]
# show_specindex(data_PWN_flag, True)
# sv_fig("../Output/SpecIndex_PWN-core_flag")

# file_PWNcore = "../DAta/i.cube.convolved.SNR.txt"
# data_PWN = np.loadtxt(file_PWNcore)
# show_specindex(data_PWN, True)
# # sv_fig("../Output/SpecIndex_PWN-core_orig")
# data_PWN_flag = data_PWN[55:]
# show_specindex(data_PWN_flag, True)
# sv_fig("../Output/SpecIndex_PWN-core_flag")

data_shell  = np.array([
    (8.354910000000E+08, 0.0425, 0.0053),
    (9.074910000000E+08, 0.0413, 0.0042),
    (9.794910000000E+08, 0.0390, 0.0036),
    (1.051490000000E+09, 0.0362, 0.0032)
])


plt.show()

# %%


reg = f'fk5;ellipse(210.1915144,-63.4279033,45.157",22.803",333.13358)'

# %%
