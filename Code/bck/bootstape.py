#%%
import numpy as np
datafile = "../Data/bootstrap.dat"

data = np.loadtxt(datafile, delimiter=",")

x = data[:, 0]  # x轴位置
y = data[:, 1]  # y轴位置
e = data[:, 2]   # y轴误差
print(f"{x, y, e}")

# fitting with polyfit

popt, cov = np.polyfit(x, y, 1, cov=True) # 1次多项式拟合，返回协方差矩阵
err = np.sqrt(np.diag(cov))  # 由协方差矩阵计算多项式系数误差
print(popt, err)         # 输出拟合结果及误差
print(np.poly1d(popt))   # 输出拟合函数

#%%
from scipy.optimize import curve_fit
f = lambda x, *p: np.polyval(p, x)   # 定义待拟合函数
# 依次传入拟合函数、拟合数据，参数初值，及误差
popt, cov = curve_fit(f, x, y, [1, 1], sigma = e, absolute_sigma = True)
err = np.sqrt(np.diag(cov))
print(popt,err)         # 输出拟合结果及误差
print(np.poly1d(popt))  # 输出拟合函数
# %%
# manual check 
import pylab as plt
xx = np.linspace(1,5,5)
aa = popt[0] ; bb = popt[1]
da = err[0] ; db = err[1]
plt.errorbar(x,y,yerr=e, fmt="o",capsize=4)  # 数据点
plt.plot(xx, np.poly1d([aa,bb])(xx),'r-')  # 绘制最佳拟合结果 
plt.plot(xx, np.poly1d([aa+da,bb+db])(xx),'c-')
plt.plot(xx, np.poly1d([aa-da,bb-db])(xx),'c-')
plt.plot(xx, np.poly1d([aa+da,bb-db])(xx),'m-')
plt.plot(xx, np.poly1d([aa-da,bb+db])(xx),'m-')
plt.show()
# %%
xi = np.linspace(1, 5, 20)   # 取样点越密，置信带越平滑
# 根据拟合结果生成1000个模拟函数
ps = np.random.multivariate_normal(popt, cov, 1000) 
# 根据模拟函数统计每个取样点上的y值

ysample = np.asarray([f(xi, *pi) for pi in ps])
print(ysample.shape)

# 在每个取样点上提取 68% 区间上下边界 
lower = np.percentile(ysample, 16, axis=0)
upper = np.percentile(ysample, 84, axis=0)
plt.errorbar(x,y,yerr=e, fmt="o",capsize=4)
plt.plot(xi, np.poly1d(popt)(xi), 'r-')  # 绘制最佳拟合函数
plt.fill_between(xi, lower, upper, color="gray",alpha=0.2)
plt.show()
# %%
import uncertainties as unc
from uncertainties import unumpy
print(cov)

a, b = unc.correlated_values(popt, cov) # 将协方差矩阵信息导入系数
print(a, b)
py = a * xi + b     # 将采样点xi带入拟合函数得到目标数据特征
print(py)
nom = unumpy.nominal_values(py)  # 得出各采样点处的期望值
std = unumpy.std_devs(py)        # 得出各采样点处的方差
plt.errorbar(x,y,yerr=e, fmt="o",capsize=4) # 绘制数据点
plt.plot(xi, nom, c='r')             #  绘制最佳拟合函数
plt.plot(xi, nom - std, c='c')       #  绘制最佳拟合的1σ下限
plt.plot(xi, nom + std, c='c')       #  绘制最佳拟合的1σ上限

# %%


# %%
