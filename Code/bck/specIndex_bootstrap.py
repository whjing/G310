
# %%
    
# ------------------------------- #
# 转换成一次多项式去拟合
# 感觉拟合出来的好好像也对？
# 还是要用定义的表达式去拟合
# 
# ------------------------------- # 
    
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties as unc
from uncertainties import unumpy

def specIndex(freq, A, alpha):
    return A * freq ** (-alpha)

def __main__():
    datafile = "../Data/freqFluxErr_SNR.dat"
    freq, flux, flux_err = np.genfromtxt(datafile, delimiter=",", unpack=True)

    freq_lg = np.log(freq)
    flux_lg = np.log(flux)
    flux_err_u = np.log10(flux + flux_err) - np.log10(flux)
    flux_err_l = abs(np.log10(flux - flux_err) - np.log10(flux))
    flux_err_lg = np.maximum(flux_err_u, flux_err_l)

    f = lambda x, *p: np.polyval(p, x)  
    popt, pcov = curve_fit(f, freq_lg, flux_lg, [1, 1], sigma = flux_err_lg, absolute_sigma = True)
    err = np.sqrt(np.diag(pcov))

    print(popt,err) 
    print(np.poly1d(popt))

    a, b = unc.correlated_values(popt, pcov)

    flux_func = a * freq_lg + b
    print(a, b)
    nom = unumpy.nominal_values(flux_func)
    std = unumpy.std_devs(flux_func)   
    plt.errorbar(freq_lg, flux_lg, yerr=flux_err_lg, fmt="o",capsize=4) # 绘制数据点
    plt.plot(freq_lg, nom, c='r',label = f"alpha = {a}")             #  绘制最佳拟合函数
    plt.plot(freq_lg, nom - std, c='c')       #  绘制最佳拟合的1σ下限
    plt.plot(freq_lg, nom + std, c='c')       #  绘制最佳拟合的1σ上限
    plt.xlabel('Frequency')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    __main__()
#%%
    
# ——————————————————————————————-
# 用公式做拟合
# —————————————————————————--————


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties as unc
from uncertainties import unumpy


# 定义幂律函数模型
def power_law(freq, A, alpha):
    return A * freq**(alpha)

def __main__():

    datafile = "../Data/freqFluxErr_SNR.dat"
    freq, flux, err = np.genfromtxt(datafile, delimiter=",", unpack=True)
    print(freq)
    f_err = 0.4
    err = err * f_err
    # 进行带有误差的幂律拟合
    popt, pcov = curve_fit(power_law, freq, flux , sigma=err, absolute_sigma = True, maxfev=10000)#[1, -0.4]
    # # 输出拟合结果
    print("拟合结果： A =", popt[0], ", alpha =", popt[1], "alpha error", pcov[1])
    A, alpha = unc.correlated_values(popt, pcov) # 将协方差矩阵信息导入系数
    print(A, alpha)
    flux_fit = power_law(freq, A, alpha)   # 将采样点xi带入拟合函数得到目标数据特征
    nom = unumpy.nominal_values(flux_fit)  # 得出各采样点处的期望值
    std = unumpy.std_devs(flux_fit)        # 得出各采样点处的方差
    plt.errorbar(freq, flux, yerr=err, fmt="o", capsize=4) # 绘制数据点
    plt.plot(freq, nom, c='r', label = f"alpha = {alpha}")     
    # plt.plot(freq, nom, c='r', label = f"alpha = {alpha}")     
    #  绘制最佳拟合函数
    plt.plot(freq, nom - std, c='c')       #  绘制最佳拟合的1σ下限
    plt.plot(freq, nom + std, c='c')       #  绘制最佳拟合的1σ上限
    # plt.loglog()
    plt.xlabel('Frequency')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()
if __name__ == "__main__":
    __main__()
#%%

    xi = np.linspace(0.8e9, 1.1e9, 20)   # 取样点越密，置信带越平滑
    # 根据拟合结果生成1000个模拟函数
    ps = np.random.multivariate_normal(popt, pcov, 10)
    for pi in ps:
        print(pi)
        ysample = np.asarray([power_law(xi, A, pi)])
        print(ysample)
        plt.plot(xi, ysample[0], 'r-', label=f"{popt[1]}")
    plt.show()
    # print(ps_alpha)
    # 根据模拟函数统计每个取样点上的y值
    # print(ps)
    # ysample = np.asarray([f(xi, *pi) for pi in ps])
    # ysample = np.asarray([power_law(xi, *pi) for pi in ps])
    # ysample = np.asarray([power_law(xi, A,alpha) for A,alpha in ps_A, ps_alpha])
    

    # 在每个取样点上提取 68% 区间上下边界 
    # lower = np.percentile(ysample, 16, axis=0)
    # upper = np.percentile(ysample, 84, axis=0)
    # plt.errorbar(freq, flux ,yerr=err, fmt="o", capsize=4)
    # plt.plot(xi, power_law(xi, popt[0], popt[1]), 'r-', label=f"{popt[1]}")  # 绘制最佳拟合函数
    # plt.fill_between(xi, lower, upper, color="gray",alpha=0.2)
    # plt.loglog()
    # plt.legend()
    # plt.show()

if __name__ == "__main__":
    __main__()
# %%
