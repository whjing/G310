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

    freq, flux, err = np.genfromtxt(datafile, delimiter=",", unpack=True) 
    freq = freq * 1e-9
    flux = flux * 1e3

    err = err * f_err * 1e3

    popt, pcov = curve_fit(power_law, freq, flux, sigma=err, absolute_sigma=True, maxfev=1000000)
    alpha_list = []
    for flux_fit in np.random.normal(flux, err, (100000, len(freq))):
        popt, pcov = curve_fit(power_law, freq, flux_fit, maxfev=10000)
        alpha_list.append(popt[1])

    alpha_array = np.array(alpha_list)
    alpha = np.round(np.mean(alpha_array), 1)
    alpha_err = np.round(np.std(alpha_array), 1)

    print(f"results: alpha = {alpha}, alpha error = {alpha_err}")
    return freq, flux, err, popt, alpha, alpha_err

def plot_specIndex(datafilePWN, datafileSNR):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    freq_lin = np.linspace(0.8,1.1,20)

    # Plot PWN data
    freq, flux, err, popt, alpha, alpha_err = bootstrap(datafilePWN, f_err=1)
    ax1.errorbar(freq, flux, yerr=err, fmt="o", capsize=4)
    ax1.plot(freq_lin, power_law(freq_lin, *popt), color='red', label=r"$\alpha = {%s} \pm {%s}$" % (alpha, alpha_err))
    ax1.legend()

    # Plot SNR data
    freq, flux, err, popt, alpha, alpha_err = bootstrap(datafileSNR, f_err=0.4)
    ax2.errorbar(freq, flux, yerr=err, fmt="o", capsize=4)
    ax2.plot(freq_lin, power_law(freq_lin, *popt), color='red', label=r"$\alpha = {%s} \pm {%s}$" % (alpha, alpha_err))
    ax2.legend()

    ax2.set_xlabel(f'frequency (GHz)')
    ax1.set_ylabel(r'flux density ($\rm{mJy}$)')
    ax2.set_ylabel(r'flux density ($\rm{mJy}$)')
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_yscale("log")

    ax1.set_xticks([0.80, 0.85, 0.9, 0.95, 1, 1.05, 1.10], [0.80, 0.85, 0.9, 0.95, 1.00, 1.05, 1.10])
    ax1.set_yticks([225, 230, 235, 240, 245, 250, 255,260], [225, 230, 235, 240, 245, 250, 255, 260]) 
    ax2.set_yticks([34, 36, 38, 40, 42, 44, 46], [34, 36, 38, 40, 42, 44, 46]) #

    # ax1.set_xticks()

    plt.tight_layout()
    sv_fig("../figures/spec2d")
    plt.show()

def main():
    datafile_PWN = "../Data/freqFluxErr_PWN.dat"
    datafile_SNR = "../Data/freqFluxErr_SNR.dat"
    plot_specIndex(datafile_PWN, datafile_SNR)

if __name__ == "__main__":
    main()

# %%
