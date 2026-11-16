"""Generate publication-style Chinese figures from archived experiment tables."""

from pathlib import Path

import matplotlib
import networkx as nx
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from run_experiments import BASE, OUT, real_graphs

FIG = BASE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
font = FIG / "NotoSansCJKsc-Regular.otf"
for system_font in font_manager.findSystemFonts():
    if "NotoSansCJK" in Path(system_font).name:
        font_manager.fontManager.addfont(system_font)
if font.exists():
    font_manager.fontManager.addfont(str(font))
    name = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams["font.family"] = name
else:
    for candidate in [
        "Noto Sans CJK SC",
        "Noto Sans CJK JP",
        "SimHei",
        "Microsoft YaHei",
    ]:
        if any(item.name == candidate for item in font_manager.fontManager.ttflist):
            plt.rcParams["font.family"] = candidate
            break
    else:
        raise RuntimeError(
            "Chinese font not found. Install Noto Sans CJK or place NotoSansCJKsc-Regular.otf in figures/."
        )
plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 190,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.unicode_minus": False,
    }
)
COL = ["#245b88", "#d37837", "#338270"]


def save(name):
    plt.tight_layout()
    plt.savefig(FIG / name, bbox_inches="tight")
    plt.close()


graphs, truth = real_graphs()
G = graphs["ca-GrQc"]
hub = max(G, key=G.degree)
H = nx.ego_graph(G, hub)
c = nx.community.louvain_communities(H, weight=None, seed=20261116)
labels = {v: i for i, x in enumerate(c) for v in x}
fig, ax = plt.subplots(figsize=(8, 5.5))
pos = nx.spring_layout(H, seed=20261116, iterations=160)
nx.draw_networkx_edges(H, pos, alpha=0.18, width=0.65, ax=ax)
nx.draw_networkx_nodes(
    H,
    pos,
    node_color=[labels[v] for v in H],
    cmap="tab20",
    node_size=[100 if v == hub else 26 for v in H],
    ax=ax,
)
ax.set_title(f"真实科研合作网络的局部结构：最高度节点的一阶邻域（{len(H)} 节点）")
ax.text(
    0,
    -0.03,
    "仅用于局部结构展示；全文统计在完整网络或其最大连通分量上计算。",
    transform=ax.transAxes,
    fontsize=9,
)
ax.axis("off")
save("q2_real_collaboration.png")

s = pd.read_csv(OUT / "q4_er_summary.csv")
df = pd.read_csv(OUT / "q4_er_replicates.csv")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.1))
for model, col in zip(["Gnp", "Gnm"], COL):
    sub = df[(df.N == 1000) & (df.p == 0.1) & (df.model == model)]
    if model == "Gnp":
        axs[0].hist(sub.m, bins=15, alpha=0.65, color=col, label="G(N,p)")
    else:
        axs[0].axvline(sub.m.iloc[0], color=col, lw=2, label="G(N,M) 固定边数")
    axs[1].scatter(sub.m, sub.avg_path, s=13, alpha=0.5, color=col, label=model)
axs[0].set(xlabel="边数 M", ylabel="重复次数", title="N=1000，p=0.1：边数波动")
axs[1].set(
    xlabel="实际边数 M",
    ylabel="最大连通分量平均距离 L",
    title="边数波动与路径长度的关系",
)
for ax in axs:
    ax.legend(frameon=False, fontsize=9)
save("q4_edge_fluctuations.png")

s = pd.read_csv(OUT / "q5_powerlaw_summary.csv")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for gamma, col in zip([2.2, 3.0, 4.0], COL):
    z = s[s.gamma == gamma]
    axs[0].errorbar(
        z.N,
        z.avg_path_mean,
        yerr=z.avg_path_sd,
        marker="o",
        capsize=3,
        color=col,
        label=f"γ={gamma:g}",
    )
    axs[1].errorbar(
        z.N,
        z.diameter_mean,
        yerr=z.diameter_sd,
        marker="s",
        capsize=3,
        color=col,
        label=f"γ={gamma:g}",
    )
for ax in axs:
    ax.set_xscale("log")
    ax.set_xlabel("节点数 N")
    ax.legend(frameon=False)
    ax.grid(alpha=0.15)
axs[0].set_ylabel("平均距离 L")
axs[1].set_ylabel("精确直径 D")
fig.suptitle("截断幂律网络：20 次独立重复的均值 ± 标准差")
save("q5_distances.png")

deg = pd.read_csv(OUT / "q5_degree_samples.csv")
fig, axs = plt.subplots(1, 3, figsize=(12, 3.7))
for ax, gamma, col in zip(axs, [2.2, 3.0, 4.0], COL):
    x = deg[(deg.N == 10000) & (deg.gamma == gamma)].degree.to_numpy()
    ks = np.arange(2, 101)
    p = ks.astype(float) ** (-gamma)
    p /= p.sum()
    ax.loglog(
        ks,
        np.array([(x >= k).mean() for k in ks]),
        ".",
        color=col,
        label="实际度序列（rep=0）",
    )
    ax.loglog(
        ks, np.cumsum(p[::-1])[::-1], "-", color="0.3", lw=1, label="目标截断分布"
    )
    ax.set(xlabel="度 k", ylabel="P(K≥k)", title=f"N=10000，γ={gamma:g}")
    ax.legend(frameon=False, fontsize=7)
save("q5_degree_validation.png")

fig, axs = plt.subplots(1, 2, figsize=(10.5, 4))
fit = pd.read_csv(OUT / "q5_empirical_tail_fit.csv")
for name, col in zip(graphs, COL):
    x = np.array([d for _, d in graphs[name].degree() if d > 0])
    ks = np.unique(x)
    axs[0].loglog(ks, [(x >= k).mean() for k in ks], ".-", ms=3, label=name, color=col)
    f = fit[fit.network == name].iloc[0]
    from scipy.special import zeta

    grid = np.arange(int(f.xmin), int(x.max()) + 1)
    axs[1].loglog(
        grid,
        zeta(f.alpha, grid.astype(float)) / zeta(f.alpha, float(f.xmin)),
        color=col,
        label=f"{name}：拟合尾部",
    )
    tail = x[x >= f.xmin]
    u = np.unique(tail)
    axs[1].loglog(u, [(tail >= k).mean() for k in u], ".", color=col, alpha=0.6)
axs[0].set(xlabel="度 k", ylabel="整体 CCDF", title="真实网络：描述性度分布")
axs[1].set(
    xlabel="度 k", ylabel="条件尾部 CCDF", title="离散幂律尾部拟合：曲线外观不代替检验"
)
for ax in axs:
    ax.legend(frameon=False, fontsize=8)
save("q5_empirical_tail.png")

s = pd.read_csv(OUT / "q6_ws_summary.csv")
fig, axs = plt.subplots(1, 2, figsize=(11, 4.1))
for ax, N in zip(axs, [100, 1000]):
    z = s[s.N == N]
    L0 = z.loc[z.p == 0, "avg_path_mean"].iloc[0]
    C0 = z.loc[z.p == 0, "clustering_mean"].iloc[0]
    ax.errorbar(
        z.p,
        z.avg_path_mean / L0,
        yerr=z.avg_path_sd / L0,
        marker="o",
        label="L(p)/L(0)",
        capsize=2,
        color=COL[0],
    )
    ax.errorbar(
        z.p,
        z.clustering_mean / C0,
        yerr=z.clustering_sd / C0,
        marker="s",
        label="C(p)/C(0)",
        capsize=2,
        color=COL[1],
    )
    ax.set_xscale("symlog", linthresh=0.0001)
    ax.set(xlabel="重连概率 p（含 p=0）", ylabel="归一化指标", title=f"N={N}，K=6")
    ax.legend(frameon=False)
    ax.grid(alpha=0.15)
save("q6_small_world.png")

K = nx.karate_club_graph()
ass = pd.read_csv(OUT / "q7_assignments.csv")
runs = pd.read_csv(OUT / "q7_community_replicates.csv")
r = (
    runs[(runs.network == "Karate") & (runs.algorithm == "Louvain")]
    .sort_values(["modularity", "rep"], ascending=[False, True])
    .iloc[0]
)
a = ass[(ass.network == "Karate") & (ass.algorithm == "Louvain") & (ass.rep == r.rep)]
pp = dict(zip(a.node, a.pred))
gt = {v: int(K.nodes[v]["club"] != "Mr. Hi") for v in K}
pos = nx.spring_layout(K, seed=20261116, weight=None)
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.4))
for ax, labs, title in [
    (axs[0], gt, "外部标签：俱乐部分裂后的两个阵营"),
    (axs[1], pp, "无权 Louvain：按最高模块度选择的一次运行"),
]:
    nx.draw_networkx(
        K,
        pos,
        node_color=[labs[v] for v in K],
        node_size=120,
        font_size=7,
        cmap="tab10",
        width=0.6,
        edge_color="0.75",
        ax=ax,
    )
    ax.set_title(title, fontsize=10)
    ax.axis("off")
save("q7_karate_comparison.png")

r = (
    runs[(runs.network == "email-Eu-core") & (runs.algorithm == "Louvain")]
    .sort_values(["modularity", "rep"], ascending=[False, True])
    .iloc[0]
)
a = ass[
    (ass.network == "email-Eu-core") & (ass.algorithm == "Louvain") & (ass.rep == r.rep)
]
ct = pd.crosstab(a.truth, a.pred)
fig, ax = plt.subplots(figsize=(9, 7))
im = ax.imshow(ct.to_numpy(), aspect="auto", cmap="Blues")
fig.colorbar(im, ax=ax, label="节点数量")
ax.set(
    xlabel="算法社区编号（无先验部门标签输入）",
    ylabel="真实部门编号",
    title="邮件网络：部门与结构社区的交叉计数",
)
ax.set_yticks(range(len(ct.index)), ct.index, fontsize=7)
save("q7_email_contingency.png")

s = pd.read_csv(OUT / "q7_resolution_summary.csv")
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
axs[0].errorbar(
    s.resolution, s.n_communities_mean, yerr=s.n_communities_sd, marker="o", capsize=3
)
axs[0].axhline(42, color="0.5", ls="--", label="42 个部门")
axs[0].legend(frameon=False)
axs[1].errorbar(
    s.resolution, s.ARI_mean, yerr=s.ARI_sd, marker="o", label="ARI", capsize=3
)
axs[1].errorbar(
    s.resolution, s.NMI_mean, yerr=s.NMI_sd, marker="s", label="NMI", capsize=3
)
axs[1].legend(frameon=False)
for ax in axs:
    ax.set_xlabel("分辨率参数 ρ")
    ax.grid(alpha=0.15)
axs[0].set_ylabel("社区数量")
axs[1].set_ylabel("与部门标签的一致性")
save("q7_resolution.png")
print("Generated 9 figures in", FIG)
