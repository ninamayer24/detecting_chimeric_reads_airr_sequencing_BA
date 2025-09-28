#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, math, random, re
from pathlib import Path
import pandas as pd
import numpy as np

DNA_RE = re.compile(r"[ACGTN]+", re.IGNORECASE)

def sanitize(seq: str) -> str:
    if not isinstance(seq, str):
        return ""
    return "".join(ch for ch in seq.upper() if ch in "ACGTN")

def pct_label(p: float) -> str:
    return f"{str(p).replace('.', '_')}p"

def choose_breakpoint_lenA(len_a: int, mid_frac: float) -> int:
    L = max(2, int(len_a))
    lo = max(1, int((1 - mid_frac) / 2 * L))
    hi = min(L - 1, int((1 + mid_frac) / 2 * L))
    if lo >= hi:
        lo, hi = 1, L - 1
    return random.randint(lo, hi)

def geometric_fit_p(sizes: np.ndarray) -> float:
    mean = float(sizes.mean())
    return 1.0 / max(mean, 1.000001)

def sample_chim_size_geometric(p_geom: float, p1: int, p2: int, lam: float) -> int:
    candidate = int(np.random.geometric(p_geom))
    max_allowed = max(1, int(math.floor(min(p1, p2) / lam)))
    return max(1, min(candidate, max_allowed))

def build_chimera_row(df, idx_a, idx_b, bp, vcol, idcol, cc_col, cc_new, chim_id_suffix, parents_meta=True):
    row_a = df.loc[idx_a].copy()
    row_b = df.loc[idx_b].copy()
    seq_a = sanitize(str(row_a.get(vcol, "")))
    seq_b = sanitize(str(row_b.get(vcol, "")))

    len_a = len(seq_a)
    len_b = len(seq_b)
    if len_a < 2 or len_b < 2:
        return None, None, None

    swapped = False
    if len_b < len_a:
        row_a, row_b = row_b, row_a
        seq_a, seq_b = seq_b, seq_a
        len_a, len_b = len_b, len_a
        swapped = True

    bp = max(1, min(bp, len_a - 1))
    chim_seq = seq_a[:bp] + seq_b[bp:len_a]

    row_a[vcol] = chim_seq
    row_a[cc_col] = int(cc_new)

    base_id = str(row_a[idcol]) if idcol in row_a else f"seq{idx_a}"
    row_a[idcol] = f"{base_id}{chim_id_suffix}"

    if parents_meta:
        orig_a_id = df.loc[idx_a][idcol] if idcol in df.columns else str(idx_a)
        orig_b_id = df.loc[idx_b][idcol] if idcol in df.columns else str(idx_b)
        row_a["_chimera"] = True
        row_a["_parentA_id"] = orig_a_id if not swapped else orig_b_id
        row_a["_parentB_id"] = orig_b_id if not swapped else orig_a_id
        row_a["_breakpoint"] = bp
        row_a["_lenA"] = len_a
        row_a["_lenB"] = len_b
        row_a["_swappedAB"] = swapped

    return row_a, len_a, swapped

def write_fasta(df: pd.DataFrame, vcol: str, idcol: str, cc_col: str, fasta_path: Path):
    with open(fasta_path, "w") as fh:
        for _, r in df.iterrows():
            sid  = str(r.get(idcol, "seq"))
            size = int(r.get(cc_col, 1))
            fh.write(f">{sid}, size={size}\n")
            seq = str(r.get(vcol, "") or "")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

def collect_input_files(paths, recursive=False):
    files = []
    for p in paths:
        pth = Path(p)
        if pth.is_dir():
            pattern = "**/*.tsv" if recursive else "*.tsv"
            files.extend(str(f) for f in pth.glob(pattern) if f.is_file())
        else:
            if pth.suffix.lower() == ".tsv" and pth.is_file():
                files.append(str(pth))
    files = sorted(set(files))
    if not files:
        raise SystemExit("[ERROR] No TSV files found.")
    return files

def main():
    ap = argparse.ArgumentParser(
        description="Inject chimeras; write FASTA with all sequences ('>sequence_id, size=consensus_count'). Sizes are drawn geometrically."
    )
    ap.add_argument("tsv", nargs="+", help="Input TSV file(s) or directories.")
    ap.add_argument("--recursive", action="store_true", help="Scan directories recursively for *.tsv")
    ap.add_argument("--percents", "-p", type=str, required=True, help="Comma-separated percent levels, e.g. 0.5,5,10")
    ap.add_argument("--outdir", "-o", default=".", help="Output root (creates simulated_<p>p subfolders).")
    ap.add_argument("--vcol", default="v_sequence_alignment", help="Sequence column.")
    ap.add_argument("--idcol", default="sequence_id", help="ID column.")
    ap.add_argument("--lambda", dest="lam", type=float, default=2.0, help="Abundance skew λ (parents ≥ λ × chimera).")
    ap.add_argument("--midfrac", type=float, default=0.6, help="Uniform breakpoint within central fraction (default 0.6 → 20–80%).")
    ap.add_argument("--seed", type=int, default=42, help="Random seed.")
    ap.add_argument("--min_parent_cc", type=int, default=1, help="Min consensus_count for parent selection.")
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    # parse percents
    percents = []
    for tok in args.percents.split(","):
        tok = tok.strip().replace("%", "")
        if tok == "":
            continue
        val = float(tok)
        if val < 0 or val > 100:
            raise SystemExit(f"[ERROR] percent out of range: {tok}")
        percents.append(val)

    inputs = collect_input_files(args.tsv, recursive=args.recursive)

    for tsv_path in inputs:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)

        if args.vcol not in df.columns:
            print(f"[WARN] {tsv_path}: missing column {args.vcol}. Skipping.")
            continue
        if args.idcol not in df.columns:
            df[args.idcol] = [f"seq{i}" for i in range(len(df))]

        # normalize consensus_count (prefer lowercase)
        canonical_cc = "consensus_count"
        if canonical_cc not in df.columns and "ConsensusCount" in df.columns:
            df[canonical_cc] = df["ConsensusCount"]
        if canonical_cc not in df.columns:
            df[canonical_cc] = 1
        df[canonical_cc] = pd.to_numeric(df[canonical_cc], errors="coerce").fillna(1).astype(int)
        args.cccol = canonical_cc

        eligible = df[df[args.cccol] >= args.min_parent_cc].index.to_list()
        if len(eligible) < 2:
            print(f"[WARN] {tsv_path}: not enough eligible parents (>=2). Skipping.")
            continue

        sizes_emp = df.loc[eligible, args.cccol].to_numpy()
        p_geom = geometric_fit_p(sizes_emp)

        base = Path(tsv_path).stem
        for pct in percents:
            label = pct_label(pct)
            out_dir = Path(args.outdir) / f"simulated_{label}"
            out_dir.mkdir(parents=True, exist_ok=True)
            out_tsv = out_dir / f"{base}_{label}.tsv"
            out_fa  = out_dir / f"{base}_{label}_all.fasta"

            n = len(df)
            n_chim = int(round(n * (pct / 100.0)))

            chim_rows = []
            tries = 0
            max_tries = max(200, max(1, n_chim) * 50)
            while len(chim_rows) < n_chim and tries < max_tries:
                tries += 1
                idx_a, idx_b = np.random.choice(eligible, size=2, replace=True)
                seq_a = sanitize(str(df.loc[idx_a, args.vcol]))
                seq_b = sanitize(str(df.loc[idx_b, args.vcol]))
                if len(seq_a) < 2 or len(seq_b) < 2:
                    continue

                bp = choose_breakpoint_lenA(len(seq_a), args.midfrac)
                p1 = int(df.loc[idx_a, args.cccol])
                p2 = int(df.loc[idx_b, args.cccol])
                cc_new = sample_chim_size_geometric(p_geom, p1, p2, args.lam)
                chim_id_suffix = f"_chim{len(chim_rows)+1}"

                row, len_a_eff, _ = build_chimera_row(
                    df, idx_a, idx_b, bp,
                    args.vcol, args.idcol, args.cccol,
                    cc_new, chim_id_suffix, parents_meta=True
                )
                if row is None:
                    continue
                if len(str(row[args.vcol])) != int(len_a_eff):
                    continue
                chim_rows.append(row)

            if len(chim_rows) < n_chim:
                print(f"[WARN] {tsv_path}: created {len(chim_rows)}/{n_chim} chimeras (constraints).")

            df_chim = pd.DataFrame(chim_rows) if chim_rows else pd.DataFrame(columns=df.columns)
            df_out  = pd.concat([df, df_chim], ignore_index=True)

            if "_chimera" not in df_out.columns:
                df_out["_chimera"] = df_out.get("_chimera", False).fillna(False)

            df_out.to_csv(out_tsv, sep="\t", index=False)
            write_fasta(df_out, args.vcol, args.idcol, args.cccol, out_fa)

            print(f"[OK] {tsv_path}: injected {len(df_chim)} ({pct}%) → {out_tsv.name}, {out_fa.name} in {out_dir}")

if __name__ == "__main__":
    main()
