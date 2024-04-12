User
# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import FuncFormatter
from globalSetting import *
set_plot_settings()

# 定义幂律函数模型
def power_law(freq, A, alpha):
    return A * freq**alpha

def bootstrap(datafile, f_err=0.4):
    print(f"Now we are fitting file: {datafile}")

    freq1, flux1= np.genfromtxt(datafile, delimiter=",", usecols=(0, 4),unpack=True) 

    # dataASKAPFile = "../Data/freqFluxErr_PWN.dat"
    # freq2, flux2 = np.genfromtxt(dataASKAPFile, delimiter=",", usecols=(0, 1),unpack=True) 
    # freq = np.hstack((freq1* 1e6, freq2))
    # print(freq)
    # flux = np.hstack((flux1, flux2))
    freq = freq1 * 1e6
    flux = flux1

    print(freq,flux)
    err = flux * 0.2

    # popt, pcov = curve_fit(power_law, freq, flux, sigma=err, absolute_sigma=True, maxfev=1000000)
    # print(popt[0])
    alpha_list = []
    for flux_fit in np.random.normal(flux, err, (100000, len(freq))):
        popt, pcov = curve_fit(power_law, freq, flux_fit, p0=[829, 0.4], maxfev=10000)
        alpha_list.append(popt[1])

    alpha_array = np.array(alpha_list)
    alpha = np.round(np.mean(alpha_array), 1)
    alpha_err = np.round(np.std(alpha_array), 1)


    print(f"results: alpha = {alpha}, alpha error = {alpha_err}")
    return freq, flux, err, popt, alpha, alpha_err

def plot_specIndex(datafile):
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 5), sharex=True)


    # Plot PWN data
    freq, flux, err, popt, alpha, alpha_err = bootstrap(datafile, f_err=1)
    freq_lin = np.linspace(np.min(freq),np.max(freq),20)
    ax1.errorbar(freq, flux, yerr=err, fmt="o", capsize=4)
    ax1.plot(freq_lin, power_law(freq_lin, *popt), color='red', label=r"$\alpha = {%s} \pm {%s}$" % (alpha, alpha_err))
    ax1.legend()

    ax1.set_xlabel(f'frequency (GHz)')
    ax1.set_ylabel(r'flux density ($\rm{mJy~beam^{-1}}$)')
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    # ax1.set_xticks([0.80, 0.85, 0.9, 0.95, 1, 1.05, 1.10], [0.80, 0.85, 0.9, 0.95, 1.00, 1.05, 1.10])
    # ax1.set_yticks([225, 230, 235, 240, 245, 250, 255,260], [225, 230, 235, 240, 245, 250, 255, 260]) 


    # ax1.set_xticks()

    plt.tight_layout()
    sv_fig("../figures/spec2D_all")
    plt.show()

def main():
    datafile= "../skyview_useful/IntFlux.data"
    plot_specIndex(datafile)

if __name__ == "__main__":
    main()

# %%
