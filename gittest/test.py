#%%
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# 假设你有频率和流量的数据
freq = [8e9, 9e9, 10e9]  # 示例频率数据
flux = [10, 20, 30]      # 示例流量数据

# 创建图形
plt.plot(freq, flux)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Flux Density')

# 自定义轴刻度标签的格式函数


# 显示图形
plt.show()

# %%
