#!/usr/bin/env python3
# summarize_chimera_all.py  → Output: sim_percent,sample,method,n_total,n_chim,n_nonchim,pct_chim
import re
import argparse
from pathlib import Path
import pandas as pd

# ---------- helpers ----------
def pct(a, b):
    try:
        a = int(a); b = int(b)
        return round(100.0 * a / b, 4) if b > 0 else 0.0
    except Exception:
        return ""

def read_text(p: Path) -> str:
    try:
        return p.read_text(errors="ignore")
    except Exception:
        return ""

def parse_percent_from_string(s: str):
    """
    Prozentangaben aus Strings ziehen, tolerant für:
      simulated_10_0p, simulated-0-1-p, 5p, 12.5%, 0,1%
    -> float oder None
    """
    if s is None:
        return None
    s = str(s).lower()

    # Muster mit 'p' am Ende (Unterstrich/Bindestrich/Punkt erlaubt)
    hits_p = re.findall(r'(\d+(?:[._-]\d+)?)\s*p(?![A-Za-z0-9])', s)
    if hits_p:
        tok = hits_p[-1].replace("_", ".").replace("-", ".")
        try:
            return float(tok)
        except ValueError:
            pass

    # echtes Prozentzeichen
    hits_pct = re.findall(r'(\d+(?:[.,]\d+)?)\s*%', s)
    if hits_pct:
        tok = hits_pct[-1].replace(",", ".")
        try:
            return float(tok)
        except ValueError:
            pass

    return None

def parse_percent_from_path(path: Path):
    """Durchsucht *alle* Pfadteile (Ordner + Dateiname) und nimmt die letzte Angabe."""
    for part in reversed(path.parts):
        val = parse_percent_from_string(part)
        if val is not None:
            return val
    return None

def percent_label(x):
    if x is None or x == "":
        return ""
    s = f"{float(x):.3f}".rstrip("0").rstrip(".")
    return s

# VSEARCH/USEARCH (Logs)
RE_TOOL_VSEARCH = re.compile(r'\bvsearch v', re.IGNORECASE)
RE_TOOL_USEARCH = re.compile(r'\busearch v', re.IGNORECASE)

RE_VS_UNIQUES = re.compile(
    r'Found\s+"?(\d+)"?\s*\(([\d\.]+)%\)\s+chimeras,\s+'
    r'(\d+)\s*\(([\d\.]+)%\)\s+non-chimeras,\s+'
    r'and\s+(\d+)\s*\(([\d\.]+)%\)\s+borderline sequences in\s+(\d+)\s+unique sequences',
    re.IGNORECASE
)
RE_VS_TOTALS = re.compile(
    r'Taking abundance information into account.*?'
    r'(\d+)\s*\(([\d\.]+)%\)\s+chimeras,\s+'
    r'(\d+)\s*\(([\d\.]+)%\)\s+non-chimeras,\s+'
    r'and\s+(\d+)\s*\(([\d\.]+)%\)\s+borderline sequences in\s+(\d+)\s+total sequences',
    re.IGNORECASE | re.DOTALL
)

RE_US_INPUT = re.compile(r'\buchime\d*_?(?:denovo|ref)\s+["\']?(.+?\.(?:fa|fasta))["\']?', re.IGNORECASE)
RE_US_LAST  = re.compile(r'100\.0%\s+Chimeras\s+(\d[\d,]*)\s+hits\s+\(([\d\.]+)%\)', re.IGNORECASE)
RE_UCHIME_CALL = re.compile(r'\buchime(\d*)_(denovo|ref)\b', re.IGNORECASE)

def sample_label_from_log(p: Path) -> str:
    return f"{p.parent.name}/{p.stem}"

def infer_variant_from_name(path: Path) -> str:
    s = path.name.lower()
    if "uchime3" in s and "denovo" in s: return "uchime_denovo3"
    if "uchime2" in s and "denovo" in s: return "uchime_denovo2"
    if "uchime_denovo" in s:             return "uchime_denovo"
    if "uchime_ref" in s:                return "uchime_ref"
    return "uchime_denovo"

def infer_variant_from_text_or_name(text: str, path: Path) -> str:
    m = RE_UCHIME_CALL.search(text)
    if m:
        num, kind = m.groups()
        if kind.lower() == "denovo":
            if num == "3": return "uchime_denovo3"
            if num == "2": return "uchime_denovo2"
            return "uchime_denovo"
        else:
            return "uchime_ref"
    return infer_variant_from_name(path)

def count_fasta_headers(fa_path: Path) -> int:
    n = 0
    try:
        with fa_path.open("r", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
    except Exception:
        pass
    return n

def collect_logs(root: Path):
    files = []
    for pat in ("*.log", "*.txt"):
        files += list(root.rglob(pat))
    return sorted([p for p in files if "uchime" in p.name.lower()])

def parse_vsearch(text: str, log_path: Path, abundance_only_percent: bool):
    rows = []
    sample = sample_label_from_log(log_path)
    variant = infer_variant_from_name(log_path)
    method = f"VSEARCH {variant}"
    sim_pct = parse_percent_from_path(log_path) or parse_percent_from_string(sample)

    m1 = RE_VS_UNIQUES.search(text)
    if m1:
        u_chim, _, u_non, _, _, _, u_tot = m1.groups()
        u_tot  = int(u_tot); u_chim = int(u_chim); u_non = int(u_non)
        rows.append({
            "sim_percent": percent_label(sim_pct),
            "sample": sample,
            "method": method,
            "n_total": u_tot,
            "n_chim": u_chim,
            "n_nonchim": u_non,
            "pct_chim": pct(u_chim, u_tot),
        })

    m2 = RE_VS_TOTALS.search(text)
    if m2:
        t_chim, _, t_non, _, _, _, t_tot = m2.groups()
        t_tot  = int(t_tot); t_chim = int(t_chim); t_non = int(t_non)
        rows.append({
            "sim_percent": percent_label(sim_pct),
            "sample": f"{sample}_Abundance",
            "method": method,
            "n_total": ("" if abundance_only_percent else t_tot),
            "n_chim":  ("" if abundance_only_percent else t_chim),
            "n_nonchim": ("" if abundance_only_percent else t_non),
            "pct_chim": pct(t_chim, t_tot),
        })
    return rows

def parse_usearch(text: str, log_path: Path):
    rows = []
    sample = sample_label_from_log(log_path)
    variant = infer_variant_from_text_or_name(text, log_path)
    method = f"USEARCH {variant}"
    sim_pct = parse_percent_from_path(log_path) or parse_percent_from_string(sample)

    n_total = ""
    m_in = RE_US_INPUT.search(text)
    if m_in:
        fa_raw = m_in.group(1).strip().strip('"\'')
        fa_path = Path(fa_raw)
        if not fa_path.is_absolute():
            fa_path = (log_path.parent / fa_path).resolve()
        n_val = count_fasta_headers(fa_path)
        n_total = n_val if n_val > 0 else ""

    matches = list(RE_US_LAST.finditer(text))
    if matches:
        m_last = matches[-1]
        chim = int(m_last.group(1).replace(",", ""))
        pct_end = float(m_last.group(2))
        non = (int(n_total) - chim) if isinstance(n_total, int) else ""

        rows.append({
            "sim_percent": percent_label(sim_pct),
            "sample": sample,
            "method": method,
            "n_total": n_total,
            "n_chim": chim,
            "n_nonchim": non,
            "pct_chim": round(pct_end, 4),
        })
    return rows

def summarize_logs(logs_root: Path, abundance_only_percent: bool):
    rows = []
    for lp in collect_logs(logs_root):
        txt = read_text(lp)
        if not txt: continue
        if RE_TOOL_VSEARCH.search(txt):
            rows.extend(parse_vsearch(txt, lp, abundance_only_percent))
        elif RE_TOOL_USEARCH.search(txt):
            rows.extend(parse_usearch(txt, lp))
    return rows

# CHMMAIRRa (TSVs) 
# Erlaube *_chim.tsv und *_chmm.tsv; IG/TR-Tags optional
RE_CHIM = re.compile(r'(?:_ig|_tr)?_(?:chim|chmm)\.tsv$', re.IGNORECASE)
RE_NONC = re.compile(r'(?:_ig|_tr)?_nonchim\.tsv$', re.IGNORECASE)

def is_chim_file(p: Path) -> bool:
    return bool(RE_CHIM.search(p.name))

def is_nonchim_file(p: Path) -> bool:
    return bool(RE_NONC.search(p.name))

def strip_suffixes_for_key(stem: str) -> str:
    # receptor-tag + chim/chmm/nonchim entfernen
    s = re.sub(r'(?:_ig|_tr)?_(?:chim|chmm)$', '', stem, flags=re.IGNORECASE)
    s = re.sub(r'(?:_ig|_tr)?_nonchim$', '', s, flags=re.IGNORECASE)
    # Prozent-Tag im Dateinamen raus (z.B. _0_1p / -0-1-p / .0.1p)
    s = re.sub(r'_[0-9]+(?:[_\.-][0-9]+)?p$', '', s, flags=re.IGNORECASE)
    return s

def count_rows(tsv_path: Path) -> int:
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str, engine="python")
        return len(df)
    except Exception:
        return 0

def summarize_chmm(chmm_root: Path):
    """
    Unterstützt Layouts:
      results_chmairra_sim*/simulated_0_1p/*.tsv   (neues Schema)
      .../*_chim.tsv & *_nonchim.tsv                (Fallback)
    """
    rows = []

    
    sim_dirs = [p for p in chmm_root.rglob("*") if p.is_dir() and p.name.lower().startswith("simulated")]
    visited = set()

    for sd in sorted(sim_dirs):
        sim_pct = parse_percent_from_string(sd.name) or parse_percent_from_path(sd)
        sim_lab = percent_label(sim_pct)

        tsvs = list(sd.glob("*.tsv"))
        groups = {}
        for t in tsvs:
            stem = t.stem
            key  = strip_suffixes_for_key(stem).lower()
            entry = groups.setdefault(key, {"chim": None, "nonchim": None, "sample_label": key})
            if is_chim_file(t):
                entry["chim"] = t
            elif is_nonchim_file(t):
                entry["nonchim"] = t

        for key, d in groups.items():
            chim_tsv = d["chim"]
            nonc_tsv = d["nonchim"]
            if chim_tsv is None and nonc_tsv is None:
                continue

            n_chim = count_rows(chim_tsv) if chim_tsv else 0
            n_non  = count_rows(nonc_tsv) if nonc_tsv else 0
            n_tot  = n_chim + n_non

            sample_label = f"{sd.name}/{d['sample_label']}"
            rows.append({
                "sim_percent": sim_lab,
                "sample": sample_label,
                "method": "CHMMAIRRa",
                "n_total": n_tot,
                "n_chim": n_chim,
                "n_nonchim": n_non,
                "pct_chim": pct(n_chim, n_tot),
            })
            if chim_tsv: visited.add(chim_tsv)
            if nonc_tsv: visited.add(nonc_tsv)

    # 2) Fallback: ohne 'simulated_*'
    chim_files  = [p for p in chmm_root.rglob("*_chim.tsv") if p not in visited]
    chim_files += [p for p in chmm_root.rglob("*_chmm.tsv") if p not in visited]
    nonc_files  = [p for p in chmm_root.rglob("*_nonchim.tsv") if p not in visited]

    chim_index = {}
    for cf in chim_files:
        key = strip_suffixes_for_key(cf.stem).lower()
        chim_index.setdefault(key, []).append(cf)

    for nf in nonc_files:
        key = strip_suffixes_for_key(nf.stem).lower()
        cands = chim_index.get(key, [])
        cf = None
        if cands:
            for c in cands:
                if c.parent == nf.parent:
                    cf = c; break
            if cf is None:
                cf = cands[0]

        n_chim = count_rows(cf) if cf else 0
        n_non  = count_rows(nf) if nf else 0
        n_tot  = n_chim + n_non

        sim_pct = parse_percent_from_path(nf) or parse_percent_from_path(cf) or parse_percent_from_string(nf.parent.name)
        sim_lab = percent_label(sim_pct)

        sample_label = f"{nf.parent.name}/{key}"
        rows.append({
            "sim_percent": sim_lab,
            "sample": sample_label,
            "method": "CHMMAIRRa",
            "n_total": n_tot,
            "n_chim": n_chim,
            "n_nonchim": n_non,
            "pct_chim": pct(n_chim, n_tot),
        })

    return rows

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(
        description="VSEARCH/USEARCH (uchime-Logs) + CHMMAIRRa (TSVs) → Tabelle: sim_percent,sample,method,n_total,n_chim,n_nonchim,pct_chim"
    )
    ap.add_argument("--logs-root", required=True, help="Ordner mit uchime-Logs (rekursiv).")
    ap.add_argument("--chmm-root", required=True, help="Ordner mit CHMMAIRRa-TSVs (rekursiv).")
    ap.add_argument("--out", default="chimera_summary_all.csv", help="Ausgabedatei (CSV).")
    ap.add_argument("--abundance-only-percent", action="store_true",
                    help="Bei VSEARCH-Abundance nur Prozent ausgeben (Zähler leer).")
    args = ap.parse_args()

    rows = []
    rows += summarize_logs(Path(args.logs_root), args.abundance_only_percent)
    rows += summarize_chmm(Path(args.chmm_root))

    cols = ["sim_percent","sample","method","n_total","n_chim","n_nonchim","pct_chim"]
    df = pd.DataFrame(rows, columns=cols)

    # sim_percent sauber sortieren
    df["_sim_num"] = pd.to_numeric(df["sim_percent"].astype(str).str.replace(",", ".", regex=False), errors="coerce")
    df = df.sort_values(["_sim_num","sample","method"]).drop(columns=["_sim_num"]).reset_index(drop=True)

    out_csv = Path(args.out).expanduser()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(df.to_string(index=False))
    print(f"\n[OK] gespeichert: {out_csv}")

if __name__ == "__main__":
    main()
