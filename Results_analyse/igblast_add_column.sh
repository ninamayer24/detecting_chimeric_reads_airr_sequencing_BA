#!/usr/bin/env bash
set -euo pipefail

BASE="/home-link/zxozk31/yard/apps/cellranger-9.0.1/vdj_IMGT_human/work/b3/8d8540db70715791dc4f812f0606a1/igblast_base"
DB="$BASE/database"
AUX="$BASE/optional_file/human_gl.aux"

export IGDATA="${IGDATA:-$HOME/ncbi-igblast-1.22.0}"
export PATH="$IGDATA:$IGDATA/bin:$PATH"

SEARCH_ROOT="/home-link/zxozk31"
PATTERN='*/vdj_annotation/06-annotate-metadata/PBMC_*_?CR_L00?_meta-pass.tsv'
OUTPUT_BASE="/home-link/zxozk31/Analyse_cons_count/igblast_out2"
mkdir -p "$OUTPUT_BASE"

tsv_to_fasta_sequence_only() {
  local tsv="$1" fasta_out="$2" sample="$3"
  awk -v samp="$sample" -F'\t' '
    BEGIN{sid=seq=vseq=0;n=0}
    NR==1{
      for(i=1;i<=NF;i++){if($i=="sequence_id")sid=i; if($i=="sequence")seq=i; if($i=="v_sequence_alignment")vseq=i}
      if(!seq && !vseq){exit 2}
      next
    }
    {
      s=(seq?$(seq):$(vseq)); gsub(/[.\-[:space:]]/,"",s); if(s=="")next
      id=(sid?$(sid):sprintf("%s_%06d",samp,++n)); gsub(/[[:space:]]/,"_",id)
      printf(">%s\n%s\n",id,s)
    }
  ' "$tsv" > "$fasta_out" || return 1
}

detect_receptor() {
  local inp="$1" ext="${inp##*.}" lowname
  lowname="$(basename "$inp" | tr '[:upper:]' '[:lower:]')"
  if [[ "$ext" == "tsv" ]]; then
    head -n 2000 "$inp" | grep -Eiq '\bTR[ABDG]|TCR' && { echo TR; return; }
    head -n 2000 "$inp" | grep -Eiq '\bIG[HKL]?|IGHV|IGKV|IGLV|BCR' && { echo IG; return; }
  else
    head -n 2000 "$inp" | grep -Eiq '^>.*(TR[ABDG]|TCR)' && { echo TR; return; }
    head -n 2000 "$inp" | grep -Eiq '^>.*(IG[HKL]?|IGH|IGK|IGL)' && { echo IG; return; }
  fi
  [[ "$lowname" =~ (tcr|tr[abdg]?|_tr_|-tr-) ]] && echo TR || echo IG
}

run_igblast_airr19() {
  local fasta="$1" receptor="$2" out_tsv="$3" seqtype dbV dbD dbJ
  if [[ "$receptor" == "TR" ]]; then
    seqtype="TCR"; dbV="$DB/imgt_human_tr_v"; dbD="$DB/imgt_human_tr_d"; dbJ="$DB/imgt_human_tr_j"
  else
    seqtype="Ig";  dbV="$DB/imgt_human_ig_v"; dbD="$DB/imgt_human_ig_d"; dbJ="$DB/imgt_human_ig_j"
  fi
  igblastn -germline_db_V "$dbV" -germline_db_D "$dbD" -germline_db_J "$dbJ" \
           -organism human -ig_seqtype "$seqtype" -domain_system imgt \
           -auxiliary_data "$AUX" -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 \
           -show_translation -num_threads 4 -query "$fasta" -outfmt 19 > "$out_tsv"
}

mapfile -d '' INPUTS < <(find "$SEARCH_ROOT" -type f -path "$PATTERN" -print0 | sort -z)
(( ${#INPUTS[@]} )) || { echo "keine Dateien gefunden: $PATTERN" >&2; exit 1; }

for INP in "${INPUTS[@]}"; do
  ext="${INP##*.}"
  base="$(basename "$INP")"
  sample="${base%.*}"
  sample_safe="${sample// /_}"
  receptor="$(detect_receptor "$INP")"
  receptor_tag="$(tr '[:upper:]' '[:lower:]' <<< "$receptor")"
  outdir="$OUTPUT_BASE/$sample_safe"
  mkdir -p "$outdir"

  if [[ "$ext" == "tsv" ]]; then
    fasta_tmp="$outdir/${sample_safe}.from_tsv.fasta"
    tsv_to_fasta_sequence_only "$INP" "$fasta_tmp" "$sample_safe" || { echo "skip: $INP"; continue; }
    [[ -s "$fasta_tmp" ]] || { echo "leer: $INP"; continue; }
    out_airr="$outdir/${sample_safe}_igblast_${receptor_tag}.tsv"
    run_igblast_airr19 "$fasta_tmp" "$receptor" "$out_airr"
  elif [[ "$ext" =~ ^(fa|fna|fasta)$ ]]; then
    out_airr="$outdir/${sample_safe}_igblast_${receptor_tag}.tsv"
    run_igblast_airr19 "$INP" "$receptor" "$out_airr"
  else
    echo "skip (ext): $INP"
  fi
done

echo "done → $OUTPUT_BASE"
