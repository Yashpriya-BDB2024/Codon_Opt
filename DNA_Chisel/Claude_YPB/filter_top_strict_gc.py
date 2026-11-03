"""
filter_top_strict_gc.py  -  Option 1: strict GC window, then best CAI.

Uses ONLY columns already in the grid output (Results sheet of human_grid.xlsx).
Does not touch the grid / run scripts.

Per gene (Name):
  1. Status == OK
  2. rare_codon_min_freq is set        (AvoidRareCodons was active: no codon rarer than the cut-off)
  3. remove duplicate DNA              (keep the first variant)
  4. Forbidden motifs left = 0
  5. Repeated 10-mers = 0
  6. Longest homopolymer <= 5
  7. GC% min (50nt win) >= 30 and GC% max (50nt win) <= 70
  8. GC% global between 45 and 55
  9. sort by CAI (highest first), keep the top N
Step 2 must come BEFORE step 3: it filters on how a variant was made, and identical DNA can
come from variants with different settings. Steps 4-8 are properties of the DNA itself.
No CAI cut-off and no codon-method restriction.

    python filter_top_strict_gc.py                 # reads outputs/human_grid.xlsx
    python filter_top_strict_gc.py --top-n 10
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DNA = "Optimized DNA Sequence"
MAX_HOMOPOLYMER = 5
WINDOW_GC_MIN, WINDOW_GC_MAX = 30.0, 70.0
GC_LO, GC_HI = 45.0, 55.0            # strict GC window (this script)
NEEDED = ["Name", "Variant", "Status", "CAI", "GC% global", "GC% min (50nt win)",
          "GC% max (50nt win)", "Longest homopolymer", "Repeated 10-mers",
          "Forbidden motifs left", "codon_method", "rare_codon_min_freq", DNA]


def eligible(res, gc_lo, gc_hi):
    """Steps 1-7. Returns (pool, same, funnel). pool keeps the grid's variant order."""
    missing = [c for c in NEEDED if c not in res.columns]
    if missing:
        sys.exit(f"Input 'Results' sheet is missing columns: {missing}")
    genes = list(res["Name"].unique())                       # input gene order
    funnel = pd.DataFrame(index=genes)
    funnel["variants run"] = res.groupby("Name").size()

    ok = res[res["Status"] == "OK"]
    funnel["status OK"] = ok.groupby("Name").size()
    ok = ok[ok["rare_codon_min_freq"].notna()]               # step 2, before removing duplicates
    funnel["rare-codon constraint on"] = ok.groupby("Name").size()
    same = (ok.groupby(["Name", DNA], sort=False)["Variant"]
            .agg(n_same="size", same_variants=lambda v: ", ".join(v.head(8))).reset_index())
    pool = ok.drop_duplicates(["Name", DNA]).reset_index(drop=True)
    funnel["distinct DNA"] = pool.groupby("Name").size()

    def apply(label, mask):
        nonlocal pool
        pool = pool[mask].copy()
        funnel[label] = pool.groupby("Name").size()

    apply("no forbidden motifs", pool["Forbidden motifs left"] == 0)
    apply("no repeated 10-mers", pool["Repeated 10-mers"] == 0)
    apply(f"homopolymer <= {MAX_HOMOPOLYMER}", pool["Longest homopolymer"] <= MAX_HOMOPOLYMER)
    apply(f"window GC {WINDOW_GC_MIN:.0f}-{WINDOW_GC_MAX:.0f}%",
          (pool["GC% min (50nt win)"] >= WINDOW_GC_MIN) & (pool["GC% max (50nt win)"] <= WINDOW_GC_MAX))
    apply(f"GC global {gc_lo:.0f}-{gc_hi:.0f}%", pool["GC% global"].between(gc_lo, gc_hi))
    pool["_g"] = pool["Name"].map({g: i for i, g in enumerate(genes)})
    return pool, same, funnel


def save_and_report(top, same, funnel, res_cols, out, top_n, extra=()):
    """top: chosen rows in final order with a 'Rank' column. Writes Excel + prints summary."""
    funnel = funnel.copy()
    funnel["selected"] = top.groupby("Name").size()
    funnel = funnel.fillna(0).astype(int)
    top = top.merge(same, on=["Name", DNA], how="left")      # keeps top's order
    first = ["Name", "Rank", "n_same", "same_variants"] + list(extra)
    rest = [c for c in res_cols if c not in ("Name", DNA)]
    top = (top[first + rest + [DNA]]
           .rename(columns={"n_same": "Variants giving same DNA",
                            "same_variants": "Same-DNA variant names (first 8)"}))
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    top.to_excel(out, sheet_name="Top", index=False)

    print(funnel.rename_axis("Name").reset_index().to_string(index=False))
    for g, r in funnel[funnel["selected"] < top_n].iterrows():
        print(f"WARNING: {g}: only {r['selected']} sequence(s) passed all gates (wanted {top_n})")
    summ = (top.assign(method=top["codon_method"].fillna("none"))
            .groupby("Name", sort=False)
            .agg(CAI_min=("CAI", "min"), CAI_max=("CAI", "max"),
                 GC_min=("GC% global", "min"), GC_max=("GC% global", "max"),
                 methods=("method", lambda s: ", ".join(f"{k}:{v}" for k, v in s.value_counts().items()))))
    print("\nSelected sequences per gene:")
    print(summ.round(3).to_string())
    print(f"\nSaved {len(top)} sequences to {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(HERE / "outputs" / "human_grid.xlsx"))
    ap.add_argument("--output", default=None, help="default: <input name>_strictGC_top<N>.xlsx")
    ap.add_argument("--top-n", type=int, default=10)
    args = ap.parse_args()

    res = pd.read_excel(args.input, sheet_name="Results")
    pool, same, funnel = eligible(res, GC_LO, GC_HI)

    # 9. best CAI first (stable sort: equal CAI keeps the grid's variant order)
    pool = pool.sort_values(["_g", "CAI"], ascending=[True, False], kind="mergesort")
    top = pool.groupby("Name", sort=False).head(args.top_n).copy()
    top["Rank"] = top.groupby("Name", sort=False).cumcount() + 1

    out = args.output or HERE / "outputs" / (
        Path(args.input).stem.replace("_grid", "") + f"_strictGC_top{args.top_n}.xlsx")
    save_and_report(top, same, funnel, list(res.columns), out, args.top_n)


if __name__ == "__main__":
    main()
