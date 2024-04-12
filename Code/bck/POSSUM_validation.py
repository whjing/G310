#%%
from bs4 import BeautifulSoup
import csv

with open("../Data/POSSUM_validation.html", "r", encoding="utf-8") as file:
    html_content = file.read()

soup = BeautifulSoup(html_content, "html.parser")
tables = soup.find_all("table")
table_sp = tables[-1]

thead = table_sp.find("thead")
# print(thead)

# 找到<tbody>标签
tbody = table_sp.find("tbody")

rows = tbody.find_all("tr")

#提取特定列的值并存储为CSV文件
with open("../Data/POSSUM_validation.reg", "w", encoding="utf-8") as f:
    # writer = csv.writer(csvfile)
    # writer.writerow(["ra_deg_cont", "dec_deg_cont"])  # 写入CSV文件的标题行
    # writer.writerow(["ra_deg_cont", "dec_deg_cont", "source", "fd_peak_fit"])
    f.write("# Region file format: DS9 version 4.1 \nglobal color=magenta dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n")
    for row in rows:
        # 提取特定列的值
        cells = row.find_all("td")
        ra_deg_cont = cells[1].get_text()  # 第2列的值
        dec_deg_cont = cells[2].get_text()  # 第3列的值
        source = cells[3].get_text()
        fd_peak_fit = cells[23].get_text()
        # 写入CSV文件
        # writer.writerow([ra_deg_cont, dec_deg_cont])
        f.write(f"circle({ra_deg_cont},{dec_deg_cont},0.01) # text={{{source}, {fd_peak_fit}}}\n")
        # print([ra_deg_cont, dec_deg_cont, source, fd_peak_fit])
    
# %%
