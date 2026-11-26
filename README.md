# 网络科学课程实验

本项目是网络科学课程作业的公开代码与实验数据仓库：<https://github.com/ChinaYi/network-science-homework>。报告正文仅列算法伪代码，完整 Python、C++ 实现、逐次实验数据及汇总结果在这里提供。

[RESULTS.md](RESULTS.md) 按报告第 2、4 至 7 题整理结果、方法限制与九张结果图。报告另行排版打印，课程课件、课件截图、原始作业文件和个人信息不放入本仓库。

## 与报告的对应关系

| 报告内容 | 完整实现 | 实验数据与结果 |
|---|---|---|
| 第 1、3 题：网络科学基础及与人工智能的关系 | 理论讨论，无额外训练实验 | 报告正文 |
| 第 2 题：科研合作网络 | `scripts/run_experiments.py` 的 `real_graphs`、`metrics` | `results/real_network_metrics.csv`、`figures/q2_real_collaboration.png` |
| 第 4 题：两种 ER 模型 | `scripts/run_experiments.py` 的 `er` | `results/q4_er_replicates.csv`、`results/q4_er_summary.csv` |
| 第 5 题：截断幂律图与真实图 | `power_graph`、`power`、`tail`、`real` | `results/q5_*`、`results/real_degree_null_*` |
| 第 6 题：WS 与小世界比较 | `ws`、`real` | `results/q6_ws_*`、与第 5 题共用的保度对照 |
| 第 7 题：社团识别与分辨率扫描 | `community`、`comm_metrics` | `results/q7_*` |

上述函数均位于 `scripts/run_experiments.py`。平均最短路与直径由 `scripts/path_metrics.cpp` 精确计算。结果表中的均值、标准差和置信区间来自已保存的 CSV，九张图与报告采用相同版本。

## 环境与运行

实际运行环境为 Python 3.11，依赖版本见 `requirements.txt`。精确最短路程序需要支持 C++17、OpenMP 的 `g++`，使用 4 个线程。已提交的 CSV 和 PNG 可以直接查看，不必重新计算。

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/verify_results.py
python network_science_homework.py --section validate
```

前一条命令核对已保存结果的样本数、汇总、度数和边数约束、节点划分覆盖范围。后一条命令编译 C++ 程序，将五类图的距离和直径与 NetworkX 精确结果比较，并检查六组幂律图的度序列约束。这两项检查不下载数据。

完整复现会重新生成结果，包含多次全源 BFS 和真实网络随机化，运行时间明显长于上述检查：

```bash
python network_science_homework.py
```

首次运行从 Stanford SNAP 下载三个公开文件，核对 SHA-256 后缓存于 `datasets/raw/`。原始数据和清洗图不提交，下载链接、校验值及预处理规则见 [DATA_SOURCES.md](DATA_SOURCES.md)。下载失败时可以把原文件放入该目录后重试，缓存仍会校验。

也可以分节执行：

```bash
python network_science_homework.py --section download
python network_science_homework.py --section er
python network_science_homework.py --section power
python network_science_homework.py --section ws
python network_science_homework.py --section real
python network_science_homework.py --section community
python network_science_homework.py --section tail
python network_science_homework.py --figures-only
```

图中文字使用 Noto Sans CJK SC。系统已安装该字体时可直接绘图；也可从 [Noto CJK 官方仓库](https://github.com/notofonts/noto-cjk) 下载 `NotoSansCJKsc-Regular.otf` 放到 `figures/`。字体许可见其仓库的 SIL Open Font License。已保存 PNG 已嵌入文字，不需要安装字体才能阅读。

## 实验设置

| 内容 | 设置 | 保存的随机重复数 |
|---|---|---:|
| ER 两种模型 | N=100、1000，p=0.01、0.1；G(N,p) 与固定边数 G(N,M) | 每组 100，共 800 |
| 截断离散幂律图 | N=100、1000、10000；指数 2.2、3、4；最小度 2、最大度 floor(sqrt(N)) | 每组 20，共 180 |
| WS 小世界 | N=100、1000，K=6，9 个重连概率 | 每组 40，共 720 |
| 真实图保度对照 | 两个真实图的最大连通分量，保持连通性与每个节点的度 | 每图 20，共 40 |
| 社团识别 | Karate 与邮件最大分量；Louvain、贪心模块度、异步标签传播；Karate 另做 GN | 随机算法每图 30 |
| 分辨率扫描 | 邮件图 Louvain，5 个分辨率 | 每组 10，共 50 |
| 真实度分布尾部 | 离散幂律 MLE、选择阈值、半参数自举并重新拟合 | 每图 199 |

主种子为 `20261116`，按实验节、参数和重复编号经 NumPy `SeedSequence` 派生。每次模拟使用的实际种子写在 CSV 中。复现应使用锁定的依赖版本；不同库版本可能改变同一种子对应的图或分区。

固定边数取 `M=floor(p*N*(N-1)/2+0.5)`，所以 N=100、p=0.01 时取 50，与独立连边模型的期望边数 49.5 相差 0.5。幂律图先拒绝不满足偶数和或可图性条件的度序列，再用 Havel–Hakimi 建图及 `20*m` 次成功双边交换，严格保持采样度数。真实图对照使用 `20*m` 次连通保度交换尝试，实际接受数另存。

所有图均按无向简单图分析。平均最短路与直径在最大连通分量上用全源 BFS 精确计算；同时记录整图连通性、孤立点和最大分量占比。最小度、最大度及平均局部聚类系数按完整分析图计算；真实图与保度对照的聚类比较则统一在对应最大分量上进行。社团算法、模块度统一忽略权重，真实标签只参与 ARI、NMI 评价。

## 文件

| 文件 | 内容 |
|---|---|
| `network_science_homework.py` | 分节与完整运行入口 |
| `scripts/run_experiments.py` | 数据预处理、生成模型、指标、社团与自举 |
| `scripts/path_metrics.cpp` | 精确全源 BFS 的 OpenMP 实现 |
| `scripts/make_figures.py` | 从结果表生成中文图 |
| `scripts/verify_results.py` | 已存结果的一致性检查 |
| `results/*_replicates.csv` | 参数、种子与每次实验指标 |
| `results/*_summary.csv` | 均值、标准差、近似 95% 置信区间半宽 |
| `results/environment.json` | 原始计算的软件版本、种子与分节运行时间 |
| `results/q7_assignments.csv` | 每次社团划分的节点标签与真实标签 |
| `results/q5_degree_samples.csv` | 各幂律参数第 0 次样本的节点度 |
| `datasets/source_checksums.json` | 下载前使用的固定 SHA-256 校验值 |
| `datasets/manifest.json` | 本次数据来源、下载时间、字节数和完整 SHA-256 |
| `figures/*.png` | 九张计算结果图 |

主要结果与方法的解释见 [RESULTS.md](RESULTS.md)。程序运行后产生的原始文件、清洗边表、编译产物和本机环境日志已通过 `.gitignore` 排除。

## 解释范围

有限截断幂律图不是无截断的渐近模型；改变指数也改变了平均度，距离差异不能全部归因于指数。双边交换的有限次数不构成已达到均匀分布的证明。真实数据的尾部检验未支持纯幂律假设，自举也没有建模节点度之间的依赖。高模块度不等于恢复真实部门，相关评价同时报告 ARI、NMI 和社区数。
