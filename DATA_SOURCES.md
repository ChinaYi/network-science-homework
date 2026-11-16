# 数据来源与预处理

## ca-GrQc 科研合作网络

- 发布者：Stanford SNAP。
- 数据说明：https://snap.stanford.edu/data/ca-GrQc.html
- 原始下载：https://snap.stanford.edu/data/ca-GrQc.txt.gz
- 时间范围：1993 年 1 月至 2003 年 4 月。
- 语义：节点是作者，共同署名形成无向合作关系。
- 清洗：合并重复和反向记录，删除 12 个自环，保留边表中出现过的全部节点。清洗后 5,242 节点、14,484 边，最大分量为 4,158 节点、13,422 边。最大分量不替代完整图，指标分别注明范围。

该数据条目引用 Leskovec, Kleinberg, Faloutsos, *Graph Evolution: Densification and Shrinking Diameters*, ACM TKDD, 2007。

## email-Eu-core 邮件网络与部门标签

- 发布者：Stanford SNAP。
- 数据说明：https://snap.stanford.edu/data/email-Eu-core.html
- 有向边表：https://snap.stanford.edu/data/email-Eu-core.txt.gz
- 部门标签：https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz
- 语义：匿名机构成员之间的邮件联系；标签包含 42 个部门。
- 清洗：按标签补全 1,005 个节点；任一方向存在邮件即建立一条无向边；合并互惠边，删除自环。清洗后 16,064 条边、19 个孤立点。社团识别使用 986 节点的最大分量，该分量仍包含全部 42 种部门标签。
- 标签只用于外部评价，不传入社团识别算法。

该数据条目引用 Yin, Benson, Leskovec, Gleich, *Local Higher-order Graph Clustering*, KDD, 2017。

## Karate Club

- 通过 `networkx.karate_club_graph()` 取得 34 节点、78 边及两个真实阵营标签。
- 文档：https://networkx.org/documentation/stable/reference/generated/networkx.generators.social.karate_club_graph.html
- 原始研究：Zachary, W. W., *An Information Flow Model for Conflict and Fission in Small Groups*, Journal of Anthropological Research, 33:452–473, 1977。
- 算法和模块度均使用 `weight=None`，不使用内置图的接触强度权重。

## 来源校验与获取

三个压缩文件的完整 SHA-256、原始 URL、字节数和实际下载时间保存在 `datasets/manifest.json`。固定校验值另存为 `datasets/source_checksums.json`，防止下载过程覆盖预期值。缓存与新下载文件均经过校验；不匹配时停止，不静默继续。

```bash
python network_science_homework.py --section download
```

本仓库不重新分发上述原始数据。请从发布者处取得，并遵循数据发布页面的使用和引用要求。

## 随机生成数据

ER、WS、截断幂律网络和保度对照均由本项目程序生成。参数、主种子及派生种子规则见 README 和 `scripts/run_experiments.py`，每次实际种子在结果 CSV 中。幂律生成使用离散截断分布、可图性拒绝采样和保度双边交换；不通过删除自环或重边来修改已采样的度序列。
