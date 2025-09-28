#!/usr/bin/env bash
set -euo pipefail

# ======= Nutzung =======
# ./run_vsearch_batch.sh <INPUT_DIR_FASTA> [OUTPUT_DIR] [THREADS]
# Beispiel:
# ./run_vsearch_batch.sh /path/to/fastas /path/to/out 8
#
# Erkennt FASTA-Endungen: .fa .fasta .fna

INPUT_DIR="/home-link/zxozk31/Analyse_cons_count/results_simulation6"
OUTPUT_DIR="/home-link/zxozk31/Analyse_cons_count/Results_vsearch_simulated6"
THREADS=4

VSEARCH="/home-link/zxozk31/Vsearch/vsearch-2.30.0-linux-x86_64/bin/vsearch"
USEARCH="/home-link/zxozk31/usearch_linux_x86_12.0-beta" 

REF_DB="/home-link/zxozk31/imgt_refs/IGHV_TRBV_combined.fasta"

# ======= Checks =======
if [[ -z "$INPUT_DIR" ]]; then
  echo "Usage: $0 <INPUT_DIR_FASTA> [OUTPUT_DIR] [THREADS]" >&2
  exit 1
fi
[[ -d "$INPUT_DIR" ]] || { echo "ERR: INPUT_DIR existiert nicht: $INPUT_DIR" >&2; exit 1; }
[[ -x "$VSEARCH" ]]   || { echo "ERR: vsearch nicht ausführbar: $VSEARCH" >&2; exit 1; }
[[ -f "$REF_DB" ]]    || { echo "ERR: Referenz-DB fehlt: $REF_DB" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

SUMMARY_FILE="${OUTPUT_DIR}/chimera_summary.tsv"
if [[ ! -f "$SUMMARY_FILE" ]]; then
  echo -e "Sample\tMethod\tChimeraPercent" > "$SUMMARY_FILE"
fi


extract_percent() {
  # liest aus vsearch-Logfile die (x.x%) nach "Taking abundance information into account"
  local logfile="$1"
  awk '
    /Taking abundance information into account/ {
      if (getline line) {
        if (match(line, /\(([0-9.]+)%\)[[:space:]]+chimeras,/, m)) {
          print m[1]; exit
        }
      }
    }
  ' "$logfile"
}


find "$INPUT_DIR" -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" \) -print0 \
| while IFS= read -r -d '' FASTA; do
 
  base="$(basename "$FASTA")"
  sample="${base%.*}"
  outdir="${OUTPUT_DIR}/${sample}_chimera_out"
  mkdir -p "$outdir"

  echo "==> Bearbeite: $FASTA"
  echo "    Sample: $sample"
  echo "    Output: $outdir"

  # 1) Gaps aus Sequenzen entfernen (., -), Header unverändert
  cleaned="${outdir}/cleaned.fa"
  awk 'BEGIN{OFS=""} /^>/ {print $0; next} {gsub(/[.-]/,""); print $0}' "$FASTA" > "$cleaned"

  # 2) Sicherstellen, dass size-Tags vorhanden sind
  uniques="${outdir}/uniques.fa"
  if grep -q ";size=" "$cleaned"; then
    cp -f "$cleaned" "$uniques"
    echo "    size-Tags gefunden → benutze vorhandene Abundanzen"
  else
    echo "    keine size-Tags → derep_fulllength + --sizeout"
    "$VSEARCH" --derep_fulllength "$cleaned" \
               --output "$uniques" \
               --sizeout --relabel_sha1 \
               --threads "$THREADS"
  fi

  # 3) Nach Größe sortieren
  sorted="${outdir}/uniques_sorted.fa"
  "$VSEARCH" --sortbysize "$uniques" \
             --output "$sorted" \
             --minsize 1 \
             --threads "$THREADS"

  # 4) uchime_denovo
  echo "    uchime_denovo"
  "$VSEARCH" --uchime_denovo "$sorted" \
             --uchimeout   "${outdir}/uchime_denovo.txt" \
             --uchimealns  "${outdir}/uchime_denovo_aln.txt" \
             --chimeras    "${outdir}/chimeras_denovo.fa" \
             --nonchimeras "${outdir}/nonchimeras_denovo.fa" \
             --threads "$THREADS" \
             > "${outdir}/uchime_denovo.log" 2>&1

  # 5) uchime2_denovo
  echo "    uchime2_denovo"
  "$VSEARCH" --uchime2_denovo "$sorted" \
             --uchimeout   "${outdir}/uchime2_denovo.txt" \
             --uchimealns  "${outdir}/uchime2_denovo_aln.txt" \
             --chimeras    "${outdir}/chimeras_uchime2.fa" \
             --nonchimeras "${outdir}/nonchimeras_uchime2.fa" \
             --threads "$THREADS" \
             > "${outdir}/uchime2_denovo.log" 2>&1

  # 6) uchime3_denovo
  echo "    uchime3_denovo"
  "$VSEARCH" --uchime3_denovo "$sorted" \
             --uchimeout   "${outdir}/uchime3_denovo.txt" \
             --uchimealns  "${outdir}/uchime3_denovo_aln.txt" \
             --chimeras    "${outdir}/chimeras_uchime3.fa" \
             --nonchimeras "${outdir}/nonchimeras_uchime3.fa" \
             --threads "$THREADS" \
             > "${outdir}/uchime3_denovo.log" 2>&1

  # 7) uchime_ref (gegen Germline-Referenz)
  echo "    uchime_ref"
  "$VSEARCH" --uchime_ref "$sorted" \
             --db "$REF_DB" \
             --uchimeout   "${outdir}/uchime_ref.txt" \
             --uchimealns  "${outdir}/uchime_ref_aln.txt" \
             --chimeras    "${outdir}/chimeras_ref.fa" \
             --nonchimeras "${outdir}/nonchimeras_ref.fa" \
             --threads "$THREADS" \
             > "${outdir}/uchime_ref.log" 2>&1

  # --- USEARCH uchime3_denovo ---
  if [[ -x "$USEARCH" ]]; then
    echo "    USEARCH uchime3_denovo"
    
    cp -f "$sorted" "$outdir/${sample}_uniques_sorted.fa"

    "$USEARCH" -uchime3_denovo "$outdir/${sample}_uniques_sorted.fa" \
               -uchimeout     "$outdir/${sample}_chimeras.txt" \
               -chimeras      "$outdir/${sample}_ch.fa" \
               -nonchimeras   "$outdir/${sample}_nonch.fa" \
               > "$outdir/${sample}_uchime3_denovo.log" 2>&1
  fi

  # 8) Zusammenfassung 
  for method in uchime_denovo uchime2_denovo uchime3_denovo uchime_ref; do
    logfile="${outdir}/${method}.log"
    percent="NA"
    if [[ -f "$logfile" ]]; then
      pval="$(extract_percent "$logfile" || true)"
      [[ -n "${pval:-}" ]] && percent="$pval"
    fi
    echo -e "${sample}\t${method}\t${percent}" >> "$SUMMARY_FILE"
  done

  echo "    DONE: $sample"
  echo "-------------------------------------------"
done

echo "Fertig. Zusammenfassung: $SUMMARY_FILE"
!/usr/bin/env bash
!/usr/bin/env bash
