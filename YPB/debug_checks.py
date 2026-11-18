"""
debug_checks.py -- verify drug_synergy_1.py step by step WITHOUT changing it.

It imports drug_synergy_1.py, loads YOUR real cache / PPI file / Drug_target.xlsx, and prints the
values going in and out of each step so you can see whether they are right.

Run:  python debug_checks.py        (Enter = accept [default])

Sections
  1 Files & network          5 Degree: what the script computes vs the true degree
  2 Distance-cache sanity    6 Random-set quality (are random sets degree-matched?)
  3 Inputs (genes / targets) 7 Hand re-computation of distances and separation
  4 Target clusters          8 Paper Eq.(1) direction vs the script's direction
  9 Optional: preview of the 'dist == 1' degree fix (in memory only, file untouched)
 10 Optional: mini end-to-end run with small n_random, every value printed
"""
import importlib.util
import inspect
import os
import random
import sys
import textwrap
import types

import networkx as nx
import numpy as np
import pandas as pd

WARNINGS = []


def ask(prompt, default="", must_exist=False):
    while True:
        shown = f" [{default}]" if default != "" else ""
        v = input(f"{prompt}{shown}: ").strip().strip('"').strip("'") or default
        if must_exist and not os.path.isfile(v):
            print(f"   !! File not found: '{v}'")
            continue
        return v


def head(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)


def ok(m):
    print(f"  [OK]   {m}")


def warn(m):
    print(f"  [WARN] {m}")
    WARNINGS.append(m)


def info(m):
    print(f"         {m}")


def import_module(path):
    spec = importlib.util.spec_from_file_location("drug_synergy_under_test", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["drug_synergy_under_test"] = mod
    spec.loader.exec_module(mod)
    return mod


def split_targets(text):          # same as the script's main(), plus stripping inner spaces
    return [t.strip() for t in str(text).split(",") if t.strip()]


def cluster(G, targets):          # same expansion as drug_synergy.main(): target + its direct neighbours
    s = set()
    for t in targets:
        if t in G:
            s |= {t} | set(G.neighbors(t))
    return s


def script_degree(sp, node):      # EXACTLY what get_degree_binning / get_degree_equivalents count
    return len([n for n, d in sp[node].items() if d == 0.0])


def patched_analyzer(ds, sp):
    """Copy of the analyzer whose degree rule counts distance == 1 (hops) instead of == 0.0. In memory only."""
    a = ds.NetworkProximityAnalyzer()
    a.shortest_paths = sp
    ns = dict(ds.__dict__)
    for name in ("get_degree_binning", "get_degree_equivalents"):
        src = textwrap.dedent(inspect.getsource(getattr(ds.NetworkProximityAnalyzer, name)))
        exec(src.replace("dist == 0.0", "dist == 1"), ns)
        setattr(a, name, types.MethodType(ns[name], a))
    return a


def randomness_quality(analyzer, G, label, nodes, n_sets=30):
    real = np.mean([G.degree(n) for n in nodes])
    sets = analyzer.get_random_nodes_matching_degree(set(nodes), n_sets, 100, None)
    rnd = np.mean([np.mean([G.degree(n) for n in s]) for s in sets])
    size = np.mean([len(s) for s in sets])
    line = (f"{label}: real set {len(nodes)} nodes, mean true degree {real:.1f} | "
            f"random sets: mean size {size:.1f}, mean true degree {rnd:.1f}")
    if real > 0 and not (0.7 <= rnd / real <= 1.3):
        warn(line + "   <-- random sets are NOT degree-matched")
    else:
        ok(line)


def main():
    head("SETUP")
    ds_path = ask("drug_synergy file", "drug_synergy_1.py", True)
    pkl = ask("distance cache (.pkl)", "shortest_paths_cache.pkl", True)
    ppi = ask("PPI file (.txt)", "STRING_PPI_filtered.txt", True)
    xlsx = ask("Excel (Module_genes + Drug_combination sheets)", "Drug_target.xlsx", True)
    sheet_mod = ask("  disease-gene sheet", "Module_genes")
    sheet_cmb = ask("  combination sheet", "Drug_combination")
    do_bins = ask("Run get_degree_binning on the FULL cache (slow, one pass over every row)? y/n", "n").lower() == "y"
    do_fix = ask("Preview the dist==1 degree fix in memory? y/n", "y").lower() == "y"
    n_mini = int(ask("Mini end-to-end run: n_random (0 = skip)", "100"))

    ds = import_module(ds_path)
    analyzer = ds.NetworkProximityAnalyzer()
    analyzer.load_shortest_paths(pkl)
    sp = analyzer.shortest_paths
    cache_nodes = analyzer.get_network_nodes()

    # ------------------------------------------------------------------ 1
    head("1. FILES & NETWORK")
    df = pd.read_csv(ppi, sep="\t")
    info(f"PPI columns: {list(df.columns)}  rows: {len(df)}")
    if not {"source", "target"} <= set(df.columns):
        sys.exit("PPI file needs 'source' and 'target' columns")
    G = nx.from_pandas_edgelist(df, "source", "target", edge_attr=True)
    info(f"graph built by drug_synergy.main(): {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    info(f"nodes in distance cache: {len(cache_nodes)}")
    if set(G.nodes) == cache_nodes:
        ok("cache nodes == PPI nodes (cache was built from this PPI file)")
    else:
        warn(f"cache and PPI differ: {len(set(G.nodes) - cache_nodes)} only in PPI, "
             f"{len(cache_nodes - set(G.nodes))} only in cache -> rebuild cache from this PPI file")
    loops = int((df["source"] == df["target"]).sum())
    (warn if loops else ok)(f"self-loops in PPI file: {loops}")
    key = df[["source", "target"]].apply(lambda r: tuple(sorted(map(str, r))), axis=1)
    dup = int(key.duplicated().sum())
    (warn if dup else ok)(f"duplicate undirected pairs in PPI file: {dup}")
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)
    info(f"connected components: {len(comps)}; largest {comps[0]}; next sizes {comps[1:5]}")
    main_cc = max(nx.connected_components(G), key=len)
    if "Weight" in df.columns:
        info("'Weight' column present but the cache is built from hops only, so it is not used")

    # ------------------------------------------------------------------ 2
    head("2. DISTANCE-CACHE SANITY (sampled)")
    rs = random.Random(0)
    sample = rs.sample(sorted(cache_nodes), min(30, len(cache_nodes)))
    (ok if all(sp[n].get(n) == 0 for n in sample) else warn)("distance of every sampled node to itself is 0")
    zeros_other = sum(1 for n in sample for m, d in sp[n].items() if d == 0 and m != n)
    (ok if zeros_other == 0 else warn)(f"entries with distance 0 other than the node itself: {zeros_other}")
    pairs = [rs.sample(sample, 2) for _ in range(40)]
    (ok if all(sp[a].get(b) == sp[b].get(a) for a, b in pairs) else warn)("table is symmetric on 40 sampled pairs")
    nb_ok = all(sp[n].get(m) == 1 for n in sample[:10] if n in G for m in G.neighbors(n))
    (ok if nb_ok else warn)("direct neighbours are at distance exactly 1 (unweighted hops table)")
    vals = [d for n in sample[:5] for d in sp[n].values()]
    info(f"example distance value types: {sorted({type(v).__name__ for v in vals})}  max distance seen: {max(vals)}")
    if any(isinstance(v, float) and v != int(v) for v in vals):
        warn("non-integer distances found: the cache looks WEIGHTED, not hops")
    good = True
    for a, b in pairs[:20]:
        if a in G and b in G:
            try:
                truth = nx.shortest_path_length(G, a, b)
            except nx.NetworkXNoPath:
                truth = None
            good &= (sp[a].get(b) == truth)
    (ok if good else warn)("cache == networkx shortest_path_length on 20 random pairs")

    # ------------------------------------------------------------------ 3
    head("3. INPUTS: disease genes and drug targets")
    modules = pd.read_excel(xlsx, sheet_name=sheet_mod)
    combos = pd.read_excel(xlsx, sheet_name=sheet_cmb)
    info(f"disease columns: {list(modules.columns)}")
    info(f"combination columns: {list(combos.columns)}")
    disease_sets = {}
    for dis in modules:
        g = [str(x).strip() for x in modules[dis].dropna()]
        miss = [x for x in g if x not in cache_nodes]
        off = [x for x in g if x in cache_nodes and x not in main_cc]
        disease_sets[dis] = set(g)
        (ok if not miss else warn)(f"{dis}: {len(g)} genes, {len(g) - len(miss)} in network, {len(miss)} missing")
        if miss:
            info("missing: " + ", ".join(miss))
        if off:
            info(f"in network but outside the main component (unreachable from most nodes): {', '.join(off)}")
    names = set(modules.columns)
    for dis in combos["Disease"].unique():
        if dis not in names:
            warn(f"Disease '{dis}' in combination sheet has no matching column in {sheet_mod} -> will be skipped")
    for _, r in combos.iterrows():
        for col in ("Drug_target1", "Drug_target2"):
            for t in split_targets(r[col]):
                if t not in cache_nodes:
                    warn(f"{r['Disease']}: target {t} NOT in the network")
                else:
                    info(f"{r['Disease']}: target {t}: true degree {G.degree(t)}, "
                         f"{'main component' if t in main_cc else 'OUTSIDE main component'}")

    # ------------------------------------------------------------------ 4
    head("4. TARGET CLUSTERS (target + all direct neighbours, as drug_synergy.main() builds them)")
    runs = []
    for _, r in combos.iterrows():
        if r["Disease"] not in disease_sets:
            continue
        ta, tb = split_targets(r["Drug_target1"]), split_targets(r["Drug_target2"])
        ca, cb = cluster(G, ta), cluster(G, tb)
        runs.append((r["Disease"], ta, tb, ca, cb))
        info(f"{r['Disease']}: A={ta} -> {len(ca)} nodes | B={tb} -> {len(cb)} nodes | shared nodes: {len(ca & cb)}")
        info("note: in a cluster every neighbour is 1 step from the target, so d_AA / d_BB come out ~1 by construction")

    # ------------------------------------------------------------------ 5
    head("5. DEGREE: what the script computes vs the true degree")
    check = sample[:6] + sorted({t for _, ta, tb, _, _ in runs for t in ta + tb if t in sp})
    print(f"  {'node':12s}{'true degree':>12s}{'script rule (dist==0.0)':>26s}")
    sc = []
    for n in check:
        sc.append(script_degree(sp, n))
        print(f"  {n:12s}{G.degree(n) if n in G else -1:>12d}{sc[-1]:>26d}")
    if len(set(sc)) == 1:
        warn(f"the script's degree is {sc[0]} for EVERY node tested: it only counts the node itself "
             f"(distance to itself = 0). Degree matching in get_degree_binning is therefore not working.")
    else:
        ok("script degrees vary between nodes")
    if do_bins:
        bins = analyzer.get_degree_binning()
        info(f"get_degree_binning on full cache -> {len(bins)} bin(s): "
             f"{[(lo, hi, len(v)) for lo, hi, v in bins][:5]}")
        if len(bins) == 1:
            warn("only ONE bin: every node is treated as a look-alike of every other node")

    # ------------------------------------------------------------------ 6
    head("6. RANDOM-SET QUALITY (analyzer.get_random_nodes_matching_degree, as the script uses it)")
    for dis, ta, tb, ca, cb in runs:
        randomness_quality(analyzer, G, f"{dis} drug A {ta}", ca)
        randomness_quality(analyzer, G, f"{dis} drug B {tb}", cb)
        randomness_quality(analyzer, G, f"{dis} disease genes", [g for g in disease_sets[dis] if g in cache_nodes])

    # ------------------------------------------------------------------ 7
    head("7. RE-COMPUTING DISTANCES / SEPARATION BY HAND")
    for dis, ta, tb, ca, cb in runs:
        D = {g for g in disease_sets[dis] if g in cache_nodes}
        script_d = analyzer.calculate_closest_distance(ca, D)
        md = nx.multi_source_dijkstra_path_length(G, D)          # independent: straight from the graph
        v = [md[a] for a in ca if a in md]
        indep = float(np.mean(v)) if v else float("inf")
        (ok if abs(script_d - indep) < 1e-9 else warn)(
            f"{dis} drug A: closest distance (mean over drug nodes of nearest disease node) "
            f"script={script_d:.4f}  independent={indep:.4f}")

        def nn_mean(S):
            S = list(S)
            vals = []
            for a in S:
                dd = nx.single_source_shortest_path_length(G, a)
                m = [dd[b] for b in S if b != a and b in dd]
                if m:
                    vals.append(min(m))
            return float(np.mean(vals)) if vals else 0.0

        def ab_mean(S, T):
            vals = []
            for X, Y in ((S, T), (T, S)):
                for a in X:
                    dd = nx.single_source_shortest_path_length(G, a)
                    m = [dd[b] for b in Y if b in dd]
                    if m:
                        vals.append(min(m))
            return float(np.mean(vals)) if vals else float("inf")
        if len(ca) + len(cb) <= 600:
            s_script = analyzer.calculate_separation(ca, cb)
            dAA, dBB, dAB = nn_mean(ca), nn_mean(cb), ab_mean(ca, cb)
            s_ind = dAB - (dAA + dBB) / 2
            (ok if abs(s_script - s_ind) < 1e-9 else warn)(
                f"{dis} separation: script s_AB={s_script:.4f}  independent: d_AA={dAA:.3f} d_BB={dBB:.3f} "
                f"d_AB={dAB:.3f} -> s_AB={s_ind:.4f}")
            s_t = analyzer.calculate_separation(set(ta), set(tb))
            info(f"{dis} separation using the TARGETS ONLY (paper definition, no neighbour expansion): "
                 f"s_AB={s_t:.4f}   <-- compare with {s_script:.4f} above")
            if (s_t < 0) != (s_script < 0):
                warn(f"{dis}: neighbour expansion CHANGES the overlap/separate decision "
                     f"(expanded {s_script:.3f} vs targets-only {s_t:.3f})")

    # ------------------------------------------------------------------ 8
    head("8. PAPER EQ.(1) vs THE SCRIPT'S DIRECTION")
    info("Paper Eq.(1): d(X,Y) = 1/||Y|| * sum over DISEASE nodes y of min over drug nodes x of d(x,y)")
    info("Script: calculate_closest_distance(drug, disease) = mean over DRUG nodes of min over disease nodes")
    for dis, ta, tb, ca, cb in runs:
        D = {g for g in disease_sets[dis] if g in cache_nodes}
        script_dir = analyzer.calculate_closest_distance(ca, D)
        paper_dir = analyzer.calculate_closest_distance(D, ca)       # swap arguments = paper direction
        info(f"{dis} drug A {ta}: script direction {script_dir:.4f} | paper direction {paper_dir:.4f}")

    # ------------------------------------------------------------------ 9
    if do_fix:
        head("9. PREVIEW: degree rule 'dist == 1' (in memory only; your file is NOT modified)")
        fixed = patched_analyzer(ds, sp)
        for dis, ta, tb, ca, cb in runs:
            randomness_quality(fixed, G, f"{dis} drug A {ta}", ca)
            randomness_quality(fixed, G, f"{dis} drug B {tb}", cb)
            randomness_quality(fixed, G, f"{dis} disease genes", [g for g in disease_sets[dis] if g in cache_nodes])

    # ------------------------------------------------------------------ 10
    if n_mini > 0:
        head(f"10. MINI END-TO-END RUN (n_random={n_mini}) - every returned value")
        for dis, ta, tb, ca, cb in runs:
            D = {g for g in disease_sets[dis] if g in cache_nodes}
            res = analyzer.classify_drug_combination_with_confidence(ca, cb, D, n_random=n_mini, alpha=0.05)
            print(f"\n  {dis}: {ta} + {tb}")
            for k, val in res.items():
                print(f"     {k:28s} {val}")

    head("SUMMARY")
    if WARNINGS:
        print(f"{len(WARNINGS)} warning(s):")
        for i, w in enumerate(WARNINGS, 1):
            print(f"  {i}. {w}")
    else:
        print("No warnings.")


if __name__ == "__main__":
    main()
