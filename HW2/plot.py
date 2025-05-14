import matplotlib.pyplot as plt
import numpy as np

# 問題名稱
labels = ['public1', 'public2', 'public3', 'public4', 'public5', 'public6']
x = np.arange(len(labels))  # x 軸的座標位置

# cluster 數據（原 method one）
cutsize_cluster = [4022, 1176, 29374, 252020, 178536, 191610]
runtime_cluster = [0.11, 0.13, 3.12, 13.48, 35.93, 49.57]

# without cluster 數據（原 method two）
cutsize_without = [20011, 2841, 103927, 342746, 1085707, 953246]
runtime_without = [0.09, 0.13, 0.41, 0.66, 1.96, 1.70]

width = 0.35  # 每個長條的寬度

# 長條圖: Cutsize 比較
fig1, ax1 = plt.subplots(figsize=(10,6))
rects1 = ax1.bar(x - width/2, cutsize_cluster, width, label='cluster')
rects2 = ax1.bar(x + width/2, cutsize_without, width, label='without cluster')

ax1.set_ylabel('Cutsize')
ax1.set_title('Cutsize 比較')
ax1.set_xticks(x)
ax1.set_xticklabels(labels)
ax1.legend()

def autolabel(rects, ax):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 pixels offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1, ax1)
autolabel(rects2, ax1)

fig1.tight_layout()

# 長條圖: Runtime 比較
fig2, ax2 = plt.subplots(figsize=(10,6))
rects1_r = ax2.bar(x - width/2, runtime_cluster, width, label='cluster')
rects2_r = ax2.bar(x + width/2, runtime_without, width, label='without cluster')

ax2.set_ylabel('Runtime (seconds)')
ax2.set_title('Runtime 比較')
ax2.set_xticks(x)
ax2.set_xticklabels(labels)
ax2.legend()

def autolabel_runtime(rects, ax):
    for rect in rects:
        height = rect.get_height()
        ax.annotate(f'{height:.2f}',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 pixels offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel_runtime(rects1_r, ax2)
autolabel_runtime(rects2_r, ax2)

fig2.tight_layout()

plt.show()

# 保存圖表
fig1.savefig("cutsize_comparison_cluster.png")
fig2.savefig("runtime_comparison_cluster.png")
print("Figures generated and saved as 'cutsize_comparison_cluster.png' and 'runtime_comparison_cluster.png'.")
