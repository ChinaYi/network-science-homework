"""Reproducible experiments. All distances and diameters are exact on the LCC.
Usage: python scripts/run_experiments.py [--section all|download|er|power|ws|real|community|tail]
No private datasets or credentials are required. Downloaded raw files are cached.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import requests
from scipy import optimize, special, stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

BASE = Path(__file__).resolve().parents[1]
RAW = BASE / "datasets/raw"
PROCESSED = BASE / "datasets/processed"
OUT = BASE / "results"
for p in [RAW, PROCESSED, OUT]:
    p.mkdir(exist_ok=True, parents=True)
MASTER = 20261116
SOURCES = {
    "ca-GrQc.txt.gz": (
        "https://snap.stanford.edu/data/ca-GrQc.txt.gz",
        "https://snap.stanford.edu/data/ca-GrQc.html",
    ),
    "email-Eu-core.txt.gz": (
        "https://snap.stanford.edu/data/email-Eu-core.txt.gz",
        "https://snap.stanford.edu/data/email-Eu-core.html",
    ),
    "email-Eu-core-department-labels.txt.gz": (
        "https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz",
        "https://snap.stanford.edu/data/email-Eu-core.html",
    ),
}


def seed(section, a=0, b=0, r=0):
    return int(
        np.random.SeedSequence(
            [MASTER, section, int(a), int(b), int(r)]
        ).generate_state(1)[0]
    )


def download():
    manifest = []
    expected = json.loads((BASE / "datasets/source_checksums.json").read_text())
    for name, (url, page) in SOURCES.items():
        p = RAW / name
        meta = p.with_suffix(p.suffix + ".json")
        if not p.exists():
            resp = requests.get(url, timeout=90)
            resp.raise_for_status()
            digest = hashlib.sha256(resp.content).hexdigest()
            if digest != expected[name]:
                raise ValueError(f"SHA-256 mismatch for downloaded {name}: {digest}")
            p.write_bytes(resp.content)
            meta.write_text(
                json.dumps(
                    {"downloaded_utc": datetime.now(timezone.utc).isoformat()}, indent=2
                )
            )
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        if digest != expected[name]:
            raise ValueError(f"SHA-256 mismatch for cached {name}: {digest}")
        rec = (
            json.loads(meta.read_text()) if meta.exists() else {"downloaded_utc": None}
        )
        rec.update(
            filename=name,
            url=url,
            source_page=page,
            bytes=p.stat().st_size,
            sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        )
        manifest.append(rec)
    (BASE / "datasets/manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf8"
    )


def read_pairs(name):
    with gzip.open(RAW / name, "rt") as f:
        return [
            tuple(map(int, x.split()[:2]))
            for x in f
            if x.strip() and not x.startswith("#")
        ]


def real_graphs():
    ca_pairs = read_pairs("ca-GrQc.txt.gz")
    ca = nx.Graph()
    ca.add_edges_from(ca_pairs)
    ca_loops = nx.number_of_selfloops(ca)
    ca.remove_edges_from(nx.selfloop_edges(ca))
    labels = dict(read_pairs("email-Eu-core-department-labels.txt.gz"))
    ep = read_pairs("email-Eu-core.txt.gz")
    ed = nx.DiGraph()
    ed.add_nodes_from(labels)
    ed.add_edges_from(ep)
    email = nx.Graph(ed)
    email.remove_edges_from(nx.selfloop_edges(email))
    graphs = {"ca-GrQc": ca, "email-Eu-core": email}
    (PROCESSED / "email_labels.json").write_text(
        json.dumps(labels, sort_keys=True), encoding="utf8"
    )
    desc = {
        "ca-GrQc": {
            "raw_rows": len(ca_pairs),
            "raw_self_loops": ca_loops,
            "raw_unique_undirected_edges_including_loops": len(
                set(tuple(sorted(e)) for e in ca_pairs)
            ),
            "preprocess": "无向简单图，合并反向/重复记录，删除自环，保留边表出现的全部节点",
        },
        "email-Eu-core": {
            "raw_rows": len(ep),
            "raw_nodes": len(ed),
            "raw_directed_edges": ed.number_of_edges(),
            "raw_self_loops": nx.number_of_selfloops(ed),
            "departments": len(set(labels.values())),
            "preprocess": "按标签补全1005节点；任一方向有邮件即连无向边；合并互惠边并删除自环；保留孤立节点",
        },
    }
    for name, G in graphs.items():
        nx.write_edgelist(G, PROCESSED / (name + ".edgelist"), data=False)
        (PROCESSED / (name + ".nodes.txt")).write_text(
            "\n".join(map(str, sorted(G))), encoding="utf8"
        )
        desc[name].update(
            n=len(G),
            m=G.number_of_edges(),
            components=nx.number_connected_components(G),
        )
    (PROCESSED / "preprocessing.json").write_text(
        json.dumps(desc, ensure_ascii=False, indent=2), encoding="utf8"
    )
    return graphs, labels


def exact_path(H):
    n = len(H)
    if n <= 1:
        return 0.0, 0
    binary = BASE / "scripts/path_metrics"
    src = BASE / "scripts/path_metrics.cpp"
    if not binary.exists() or binary.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(
            ["g++", "-O3", "-std=c++17", "-fopenmp", str(src), "-o", str(binary)],
            check=True,
        )
    nodes = list(H)
    idx = {v: i for i, v in enumerate(nodes)}
    text = f"{n} {H.number_of_edges()}\n" + "".join(
        f"{idx[u]} {idx[v]}\n" for u, v in H.edges()
    )
    p = subprocess.run(
        [str(binary)], input=text, text=True, capture_output=True, check=True
    )
    L, D, *_ = p.stdout.split()
    return float(L), int(D)


def metrics(G):
    H = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    L, D = exact_path(H)
    degrees = np.array([G.degree(v) for v in G], dtype=float)
    return dict(
        n=len(G),
        m=G.number_of_edges(),
        n_lcc=len(H),
        m_lcc=H.number_of_edges(),
        lcc_fraction=len(H) / len(G),
        components=nx.number_connected_components(G),
        isolates=nx.number_of_isolates(G),
        connected=int(nx.is_connected(G)),
        mean_degree=float(degrees.mean()),
        min_degree=int(degrees.min()),
        max_degree=int(degrees.max()),
        avg_path=L,
        diameter=D,
        clustering=nx.average_clustering(G),
        transitivity=nx.transitivity(G),
        path_method="exact_all_sources_BFS_on_LCC",
    )


def save_summary(df, groups, fields, name):
    ag = {}
    for f in fields:
        ag[f + "_mean"] = (f, "mean")
        ag[f + "_sd"] = (f, "std")
    ag["reps"] = (fields[0], "size")
    sm = df.groupby(groups, dropna=False).agg(**ag).reset_index()
    for f in fields:
        sm[f + "_ci95_half"] = (
            stats.t.ppf(0.975, sm.reps - 1) * sm[f + "_sd"] / np.sqrt(sm.reps)
        )
    sm.to_csv(OUT / name, index=False)
    return sm


def er(R=100):
    rows = []
    for N in [100, 1000]:
        for p in [0.01, 0.1]:
            S = N * (N - 1) // 2
            M = int(math.floor(S * p + 0.5))
            for model in ["Gnp", "Gnm"]:
                for r in range(R):
                    s = seed(4, N, int(1000 * p) + (0 if model == "Gnp" else 1), r)
                    G = (
                        nx.fast_gnp_random_graph(N, p, seed=s)
                        if model == "Gnp"
                        else nx.gnm_random_graph(N, M, seed=s)
                    )
                    z = metrics(G)
                    z.update(
                        N=N,
                        p=p,
                        model=model,
                        rep=r,
                        seed=s,
                        target_m_float=S * p,
                        fixed_m=M,
                    )
                    rows.append(z)
                print("ER", N, p, model, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q4_er_replicates.csv", index=False)
    save_summary(
        df,
        ["N", "p", "model"],
        [
            "m",
            "avg_path",
            "diameter",
            "min_degree",
            "max_degree",
            "lcc_fraction",
            "connected",
            "isolates",
        ],
        "q4_er_summary.csv",
    )


def power_graph(N, gamma, s):
    rng = np.random.default_rng(s)
    kmax = int(np.floor(np.sqrt(N)))
    ks = np.arange(2, kmax + 1)
    pr = ks.astype(float) ** (-gamma)
    pr /= pr.sum()
    attempts = 0
    while True:
        seq = rng.choice(ks, size=N, p=pr).astype(int).tolist()
        attempts += 1
        if sum(seq) % 2 == 0 and nx.is_graphical(seq):
            break
    G = nx.havel_hakimi_graph(seq)
    swaps = 20 * G.number_of_edges()
    nx.double_edge_swap(
        G, nswap=swaps, max_tries=100 * swaps, seed=int(rng.integers(2**32))
    )
    assert [G.degree(i) for i in range(N)] == seq
    assert not nx.number_of_selfloops(G)
    obs = np.array(seq)
    fit = optimize.minimize_scalar(
        lambda a: a * np.log(obs).sum() + N * np.log(np.sum(ks.astype(float) ** (-a))),
        bounds=(1.01, 8),
        method="bounded",
    ).x
    # Empirical vs specified *truncated* target CDF, not a goodness-of-fit p-value.
    D = float(np.max(np.abs(np.array([(obs <= k).mean() for k in ks]) - np.cumsum(pr))))
    return G, dict(
        kmin=2,
        kmax=kmax,
        sampling_attempts=attempts,
        swaps=swaps,
        gamma_mle=float(fit),
        degree_cdf_distance=D,
    )


def power(R=20):
    rows = []
    deg_rows = []
    for N in [100, 1000, 10000]:
        for gamma in [2.2, 3.0, 4.0]:
            for r in range(R):
                s = seed(5, N, int(gamma * 10), r)
                G, extra = power_graph(N, gamma, s)
                z = metrics(G)
                z.update(extra)
                z.update(N=N, gamma=gamma, rep=r, seed=s)
                rows.append(z)
                if r == 0:
                    nx.write_edgelist(
                        G, PROCESSED / f"power_N{N}_g{gamma}_rep0.edgelist", data=False
                    )
                    deg_rows.extend(
                        {"N": N, "gamma": gamma, "node": v, "degree": d}
                        for v, d in G.degree()
                    )
            pd.DataFrame(rows).to_csv(OUT / "q5_powerlaw_replicates.csv", index=False)
            print("POWER", N, gamma, flush=True)
    df = pd.DataFrame(rows)
    save_summary(
        df,
        ["N", "gamma"],
        [
            "avg_path",
            "diameter",
            "mean_degree",
            "max_degree",
            "clustering",
            "lcc_fraction",
            "gamma_mle",
            "degree_cdf_distance",
        ],
        "q5_powerlaw_summary.csv",
    )
    pd.DataFrame(deg_rows).to_csv(OUT / "q5_degree_samples.csv", index=False)


def ws(R=40):
    rows = []
    ps = [0, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
    for N in [100, 1000]:
        for p in ps:
            for r in range(R):
                s = seed(6, N, int(p * 10000), r)
                G = nx.watts_strogatz_graph(N, 6, p, seed=s)
                z = metrics(G)
                z.update(N=N, K=6, p=p, rep=r, seed=s)
                rows.append(z)
            print("WS", N, p, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "q6_ws_replicates.csv", index=False)
    save_summary(
        df,
        ["N", "K", "p"],
        ["avg_path", "clustering", "connected", "lcc_fraction"],
        "q6_ws_summary.csv",
    )


def real(R=20):
    graphs, _ = real_graphs()
    rows = []
    null = []
    for name, G in graphs.items():
        z = metrics(G)
        z["network"] = name
        rows.append(z)
        H = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        for r in range(R):
            J = H.copy()
            s = seed(8, 0 if name == "ca-GrQc" else 1, 0, r)
            # Preserve connectivity as well as every degree; independent chains start at H.
            attempts = 20 * H.number_of_edges()
            accepted = nx.connected_double_edge_swap(J, nswap=attempts, seed=s)
            assert nx.is_connected(J) and dict(J.degree()) == dict(H.degree())
            z = metrics(J)
            z.update(
                network=name,
                rep=r,
                seed=s,
                null="connected_degree_preserving_swaps",
                swap_attempts=attempts,
                accepted_swaps=accepted,
            )
            null.append(z)
            if r == 0:
                nx.write_edgelist(
                    J, PROCESSED / f"{name}_degree_null_rep0.edgelist", data=False
                )
        print("REAL+NULL", name, flush=True)
    pd.DataFrame(rows).to_csv(OUT / "real_network_metrics.csv", index=False)
    df = pd.DataFrame(null)
    df.to_csv(OUT / "real_degree_null_replicates.csv", index=False)
    save_summary(
        df,
        ["network"],
        ["avg_path", "diameter", "clustering", "transitivity"],
        "real_degree_null_summary.csv",
    )


def comm_metrics(G, comms, truth):
    comms = [set(x) for x in comms]
    idx = {v: i for i, c in enumerate(comms) for v in c}
    ns = sorted(G)
    assert len(idx) == len(G) and set(idx) == set(G)
    return dict(
        n_communities=len(comms),
        modularity=nx.community.modularity(G, comms, weight=None),
        ARI=adjusted_rand_score([truth[v] for v in ns], [idx[v] for v in ns]),
        NMI=normalized_mutual_info_score(
            [truth[v] for v in ns], [idx[v] for v in ns], average_method="arithmetic"
        ),
    )


def community(R=30):
    graphs, labels = real_graphs()
    K = nx.karate_club_graph()
    nx.set_edge_attributes(K, 1, "weight")
    E = graphs["email-Eu-core"]
    E = E.subgraph(max(nx.connected_components(E), key=len)).copy()
    graphs = {"Karate": K, "email-Eu-core": E}
    all_rows = []
    assignments = []
    gn_rows = []
    resolution_rows = []
    for name, G in graphs.items():
        truth = (
            {v: int(G.nodes[v]["club"] != "Mr. Hi") for v in G}
            if name == "Karate"
            else labels
        )
        # Community analysis uses the email LCC (986 nodes); labels are used ONLY for evaluation.
        for alg in ["Louvain", "Greedy", "LPA"]:
            for r in range(1 if alg == "Greedy" else R):
                s = seed(
                    7,
                    0 if name == "Karate" else 1,
                    {"Louvain": 0, "Greedy": 1, "LPA": 2}[alg],
                    r,
                )
                if alg == "Louvain":
                    comms = nx.community.louvain_communities(G, weight=None, seed=s)
                elif alg == "Greedy":
                    comms = nx.community.greedy_modularity_communities(G, weight=None)
                else:
                    comms = list(
                        nx.community.asyn_lpa_communities(G, weight=None, seed=s)
                    )
                z = comm_metrics(G, comms, truth)
                z.update(network=name, algorithm=alg, rep=r, seed=s)
                all_rows.append(z)
                assignments.extend(
                    dict(
                        network=name,
                        algorithm=alg,
                        rep=r,
                        node=v,
                        pred=i,
                        truth=truth[v],
                    )
                    for i, c in enumerate(comms)
                    for v in c
                )
        if name == "Karate":
            for comms in nx.community.girvan_newman(G):
                z = comm_metrics(G, comms, truth)
                z.update(network=name, algorithm="Girvan-Newman", rep=0, seed=-1)
                gn_rows.append(z)
                assignments.extend(
                    dict(
                        network=name,
                        algorithm="Girvan-Newman",
                        rep=len(comms),
                        node=v,
                        pred=i,
                        truth=truth[v],
                    )
                    for i, c in enumerate(comms)
                    for v in c
                )
                if len(comms) >= 10:
                    break
            best = max(gn_rows, key=lambda x: x["modularity"])
            all_rows.append(best)
        else:
            for gamma in [0.5, 1.0, 1.5, 2.0, 3.0]:
                for r in range(10):
                    s = seed(9, int(gamma * 10), 0, r)
                    c = nx.community.louvain_communities(
                        G, weight=None, resolution=gamma, seed=s
                    )
                    z = comm_metrics(G, c, truth)
                    z.update(
                        resolution=gamma,
                        rep=r,
                        seed=s,
                        objective=nx.community.modularity(
                            G, c, weight=None, resolution=gamma
                        ),
                    )
                    resolution_rows.append(z)
        print("COMMUNITY", name, flush=True)
    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / "q7_community_replicates.csv", index=False)
    # Deterministic algorithms have one observation: sd/CI remain missing, not zero uncertainty.
    save_summary(
        df,
        ["network", "algorithm"],
        ["n_communities", "modularity", "ARI", "NMI"],
        "q7_community_summary.csv",
    )
    pd.DataFrame(assignments).to_csv(OUT / "q7_assignments.csv", index=False)
    pd.DataFrame(gn_rows).to_csv(OUT / "q7_gn_hierarchy.csv", index=False)
    df = pd.DataFrame(resolution_rows)
    df.to_csv(OUT / "q7_resolution_replicates.csv", index=False)
    save_summary(
        df,
        ["resolution"],
        ["n_communities", "modularity", "ARI", "NMI"],
        "q7_resolution_summary.csv",
    )


def fit_tail(x, min_tail=100):
    x = np.asarray(x, dtype=int)
    vals = np.unique(x)
    best = None
    for xmin in vals:
        tail = x[x >= xmin]
        n = len(tail)
        if n < min_tail or len(np.unique(tail)) < 3:
            continue
        logs = np.log(tail).sum()
        opt = optimize.minimize_scalar(
            lambda a: a * logs + n * np.log(special.zeta(a, float(xmin))),
            bounds=(1.05, 8),
            method="bounded",
        )
        a = float(opt.x)
        grid = np.arange(xmin, int(tail.max()) + 1)
        hist = np.bincount(tail, minlength=int(tail.max()) + 1)
        ecdf = np.cumsum(hist)[grid] / n
        cdf = 1 - special.zeta(a, grid.astype(float) + 1) / special.zeta(a, float(xmin))
        D = float(np.max(abs(ecdf - cdf)))
        if best is None or D < best["KS"]:
            best = dict(xmin=int(xmin), alpha=a, n_tail=n, KS=D)
    return best


def sample_tail(rng, alpha, xmin, n):
    out = []
    while len(out) < n:
        # Conditional rejection from numpy's unbounded Zipf law; no artificial cutoff.
        batch = rng.zipf(alpha, size=max(500, 10 * (n - len(out))))
        out.extend(batch[batch >= xmin].tolist())
    return np.array(out[:n], dtype=int)


def tail(B=199):
    graphs, _ = real_graphs()
    records = []
    boots = []
    for name, G in graphs.items():
        x = np.array([d for _, d in G.degree() if d > 0])
        f = fit_tail(x)
        if f is None:
            continue
        rng = np.random.default_rng(seed(10, 0 if name == "ca-GrQc" else 1))
        body = x[x < f["xmin"]]
        tailprob = f["n_tail"] / len(x)
        for b in range(B):
            nt = int(rng.binomial(len(x), tailprob))
            synth = np.concatenate(
                [
                    rng.choice(body, size=len(x) - nt, replace=True),
                    sample_tail(rng, f["alpha"], f["xmin"], nt),
                ]
            )
            bf = fit_tail(synth)
            boots.append(
                dict(
                    network=name,
                    bootstrap=b,
                    KS=bf["KS"],
                    alpha=bf["alpha"],
                    xmin=bf["xmin"],
                )
            )
        exceed = sum(v["KS"] >= f["KS"] for v in boots if v["network"] == name)
        f.update(
            network=name,
            n_nonzero=len(x),
            bootstraps=B,
            p_bootstrap=(1 + exceed) / (B + 1),
            p_mc_se=math.sqrt(
                ((1 + exceed) / (B + 1)) * (1 - (1 + exceed) / (B + 1)) / (B + 1)
            ),
        )
        records.append(f)
        print("TAIL", f, flush=True)
    pd.DataFrame(records).to_csv(OUT / "q5_empirical_tail_fit.csv", index=False)
    pd.DataFrame(boots).to_csv(OUT / "q5_empirical_tail_bootstrap.csv", index=False)


def validate():
    checks = []
    for G in [
        nx.path_graph(9),
        nx.cycle_graph(12),
        nx.complete_graph(8),
        nx.karate_club_graph(),
        nx.gnp_random_graph(80, 0.15, seed=42),
    ]:
        H = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        L, D = exact_path(H)
        assert abs(L - nx.average_shortest_path_length(H)) < 1e-10 and D == nx.diameter(
            H
        )
        checks.append(
            dict(
                n=len(H),
                m=H.number_of_edges(),
                path=L,
                diameter=D,
                status="exactly_matches_networkx",
            )
        )
    for N in [100, 1000]:
        for g in [2.2, 3.0, 4.0]:
            G, _ = power_graph(N, g, seed(15, N, int(g * 10)))
            assert len(G) == N
    (OUT / "validation.json").write_text(json.dumps(checks, indent=2))
    print(
        "VALIDATED",
        len(checks),
        "exact BFS comparisons and 6 degree-preservation cases",
        flush=True,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--section",
        choices=[
            "all",
            "download",
            "er",
            "power",
            "ws",
            "real",
            "community",
            "tail",
            "validate",
        ],
        default="all",
    )
    args = ap.parse_args()
    start = time.time()
    if args.section != "validate":
        download()
    funcs = {
        "er": er,
        "power": power,
        "ws": ws,
        "real": real,
        "community": community,
        "tail": tail,
        "validate": validate,
    }
    todo = (
        list(funcs)
        if args.section == "all"
        else [args.section]
        if args.section != "download"
        else []
    )
    for k in todo:
        funcs[k]()
    env = dict(
        master_seed=MASTER,
        python=sys.version,
        platform=platform.platform(),
        networkx=nx.__version__,
        numpy=np.__version__,
        pandas=pd.__version__,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        section=args.section,
        elapsed_seconds=time.time() - start,
        path_method="exact_all_sources_BFS; 4 OpenMP threads",
    )
    (OUT / ("run_" + args.section + ".json")).write_text(json.dumps(env, indent=2))
    print("DONE", args.section, env["elapsed_seconds"], flush=True)
