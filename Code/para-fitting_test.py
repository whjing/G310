#%%
import numpy as np
import matplotlib.pyplot as plt

# 原始数据
data = np.array([
    [1, 2],
    [2, 3],
    [4, 2]
])

x = data[:, 0]
y = data[:, 1]

# 使用 numpy 的 polyfit 函数进行二次多项式拟合
coefficients = np.polyfit(x, y, 2)
a, b, c = coefficients

# 创建一组新的 x 值来绘制拟合的抛物线
x_fit = np.linspace(min(x), max(x), 100)
y_fit = a * x_fit ** 2 + b * x_fit + c

# 绘制原始数据点和拟合的抛物线
plt.scatter(x, y, label='Original Data')
plt.plot(x_fit, y_fit, color='red', label=f'Quadratic Fit: {a:.2f}x^2 + {b:.2f}x + {c:.2f}')

# 添加标题和标签
plt.title('Parabolic Fit')
plt.xlabel('x')
plt.ylabel('y')

# 添加图例
plt.legend()

# 显示图形
plt.show()

# 计算峰值对应的 x 值
x_peak = -b / (2 * a)

# 计算峰值对应的 y 值
y_peak = a * x_peak**2 + b * x_peak + c

# 输出拟合的系数和峰值对应的 x 和 y 值
print(f"Quadratic Coefficients: a = {a:.2f}, b = {b:.2f}, c = {c:.2f}")
print(f"Peak X Value: {x_peak:.2f}")
print(f"Peak Y Value: {y_peak:.2f}")


# %%
