#%%

import numpy as np
from io import StringIO

data = []

with open('../Fits_4-subband-useful/IntFlux.data', 'r') as f:
    while True:
        # 读取四行数据
        lines = [next(f) for _ in range(4)]
        if not lines:
            break
        
        # 合并四行数据为单个字符串
        combined_lines = ''.join(lines)

        # 使用StringIO将字符串转换为文件对象
        data_file = StringIO(combined_lines)

        # 从文件对象中提取后两列数据
        int_flux = np.loadtxt(data_file, delimiter=',', usecols=(-2, -1), max_rows=2)
        int_flux_bg = np.loadtxt(data_file, delimiter=',', usecols=(-2, -1), max_rows=2)

        # 将提取的数据存入列表
        data.append((int_flux, int_flux_bg))

print(data)


# %%
