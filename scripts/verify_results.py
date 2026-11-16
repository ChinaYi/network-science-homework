"""Check archived results without downloading data or rerunning large experiments."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "results"


def read(name):
    return pd.read_csv(OUT / name)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def compare_summary(stem, groups, expected_count):
    raw = read(stem + "_replicates.csv")
    summary = read(stem + "_summary.csv")
    require(len(raw) == expected_count, f"{stem}: replicate count")
    require(not raw.duplicated(groups + ["rep"]).any(), f"{stem}: duplicate run keys")
    for _, row in summary.iterrows():
        block = raw
        for group in groups:
            block = block[block[group] == row[group]]
        require(len(block) == int(row.reps), f"{stem}: group size")
        for column in summary.columns:
            if not column.endswith("_mean"):
                continue
            field = column[:-5]
            require(np.isclose(block[field].mean(), row[column]), f"{stem}: {column}")
            sd = block[field].std(ddof=1)
            require(
                np.isclose(sd, row[field + "_sd"], equal_nan=True),
                f"{stem}: {field} SD",
            )
            ci = stats.t.ppf(0.975, len(block) - 1) * sd / np.sqrt(len(block))
            require(
                np.isclose(ci, row[field + "_ci95_half"], equal_nan=True),
                f"{stem}: {field} CI",
            )
    return raw


def main():
    er = compare_summary("q4_er", ["N", "p", "model"], 800)
    power = compare_summary("q5_powerlaw", ["N", "gamma"], 180)
    ws = compare_summary("q6_ws", ["N", "K", "p"], 720)
    null = compare_summary("real_degree_null", ["network"], 40)
    community = compare_summary("q7_community", ["network", "algorithm"], 123)
    compare_summary("q7_resolution", ["resolution"], 50)

    fixed = er[er.model == "Gnm"]
    require((fixed.m == fixed.fixed_m).all(), "ER: fixed edge count")
    require((ws.m == ws.N * ws.K // 2).all(), "WS: edge count")
    require((power.min_degree >= 2).all(), "Power: minimum degree")
    require((power.max_degree <= power.kmax).all(), "Power: degree cutoff")
    require((null.connected == 1).all(), "Null graphs: connected")
    for frame in [er, power, ws, null]:
        require(np.allclose(frame.mean_degree, 2 * frame.m / frame.n), "Handshaking")
        require((frame.avg_path <= frame.diameter + 1e-12).all(), "L <= diameter")

    samples = read("q5_degree_samples.csv")
    for (n, gamma), group in samples.groupby(["N", "gamma"]):
        row = power[(power.N == n) & (power.gamma == gamma) & (power.rep == 0)].iloc[0]
        require(len(group) == n and group.node.nunique() == n, "Degree sample coverage")
        require(group.degree.sum() == 2 * row.m, "Degree sample sum")

    assignments = read("q7_assignments.csv")
    for (network, algorithm, rep), group in assignments.groupby(
        ["network", "algorithm", "rep"]
    ):
        expected = 34 if network == "Karate" else 986
        require(len(group) == expected, "Community assignment count")
        require(group.node.nunique() == expected, "Community assignment uniqueness")
    for _, row in community.iterrows():
        rep = int(row.n_communities) if row.algorithm == "Girvan-Newman" else row.rep
        group = assignments[
            (assignments.network == row.network)
            & (assignments.algorithm == row.algorithm)
            & (assignments.rep == rep)
        ]
        ari = adjusted_rand_score(group.truth, group.pred)
        nmi = normalized_mutual_info_score(
            group.truth, group.pred, average_method="arithmetic"
        )
        require(np.isclose(ari, row.ARI) and np.isclose(nmi, row.NMI), "ARI/NMI")

    bootstrap = read("q5_empirical_tail_bootstrap.csv")
    for _, row in read("q5_empirical_tail_fit.csv").iterrows():
        group = bootstrap[bootstrap.network == row.network]
        require(len(group) == row.bootstraps == 199, "Bootstrap replicate count")
        p = (1 + (group.KS >= row.KS).sum()) / (len(group) + 1)
        require(np.isclose(p, row.p_bootstrap), "Corrected bootstrap p-value")

    manifest_path = OUT / "artifact_manifest.json"
    if manifest_path.exists():
        for filename, expected in json.loads(manifest_path.read_text()).items():
            digest = hashlib.sha256((BASE / filename).read_bytes()).hexdigest()
            require(digest == expected, f"Artifact SHA-256: {filename}")
    print(
        "PASS: archived counts, summaries, constraints, partitions, bootstrap and hashes"
    )


if __name__ == "__main__":
    main()
