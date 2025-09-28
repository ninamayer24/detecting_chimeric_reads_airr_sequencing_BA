using CHMMAIRRa

println("Starte CHMMAIRRa (Batch, rekursiv)")

# Pfade
const IN_DIR     = "/home-link/zxozk31/Analyse_cons_count/results_simulation6"   
const OUT_DIR    = "/home-link/zxozk31/Analyse_cons_count/results_chmairra_sim6"  
const VFA_IGHV   = "/home-link/zxozk31/imgt_refs/IGHV_clean_simple.fasta"
const VFA_TRBV   = "/home-link/zxozk31/imgt_refs/TRBV_clean_simple.fasta"
const MAFFT_BIN  = "/home-link/zxozk31/.miniconda3/envs/chimera-env/bin/mafft"


const WRITE_NONCHIM  = true
const WRITE_CHIM     = true
const WRITE_CHIM_ALN = true

#2) Receptor & Referenz aus Dateinamen ableiten
function infer_receptor_and_vfa(path::AbstractString)
    s = lowercase(basename(path))
    if occursin("tcr", s) || occursin("trbv", s) || occursin("trb", s)
        return ("TCR", VFA_TRBV, "tr")
    else
        return ("IG",  VFA_IGHV, "ig")
    end
end

# 3) Alle .tsv rekursiv einsammeln 
function collect_tsvs(root::AbstractString)
    tsvs = String[]
    for (dir, _subdirs, files) in walkdir(root)
        for f in files
            if endswith(lowercase(f), ".tsv")
                push!(tsvs, joinpath(dir, f))
            end
        end
    end
    sort(tsvs)
end

tsv_paths = collect_tsvs(IN_DIR)

if isempty(tsv_paths)
    println("Keine .tsv gefunden in: ", IN_DIR)
else
    for tsv in tsv_paths
        receptor, v_fasta, tag = infer_receptor_and_vfa(tsv)

        # Name des direkten Elternordners (z. B. "simulated_0_5p" oder "simulated_10p")
        sim_folder = basename(dirname(tsv))

        # Zielordner: OUT_DIR/<simulated_xxx>
        outdir = joinpath(OUT_DIR, sim_folder)
        mkpath(outdir)

        base   = basename(tsv)                      # z.B. sampleX_0_5p.tsv
        prefix = replace(base, r"\.tsv$" => "")     # ohne .tsv

        out_main     = joinpath(outdir, "$(prefix)_$(tag)_chmm.tsv")
        out_nonchim  = WRITE_NONCHIM  ? joinpath(outdir, "$(prefix)_$(tag)_nonchim.tsv")           : nothing
        out_chim     = WRITE_CHIM     ? joinpath(outdir, "$(prefix)_$(tag)_chim.tsv")              : nothing
        out_chim_aln = WRITE_CHIM_ALN ? joinpath(outdir, "$(prefix)_$(tag)_chim_alignments.fasta") : nothing

        println("-> ", sim_folder, " | ", basename(tsv), "  (", receptor, ")")
        println("   Assignments: ", out_main)
        if WRITE_NONCHIM  println("   Non-chimeric-MiAIRR: ", out_nonchim) end
        if WRITE_CHIM     println("   Chimeric-MiAIRR:     ", out_chim) end
        if WRITE_CHIM_ALN println("   Chim. Alignments:    ", out_chim_aln) end

        detect_chimeras_from_files(
            v_fasta,
            tsv,
            out_main;
            receptor = receptor,
            mafft    = MAFFT_BIN,
            non_chimeric_MiAIRR = out_nonchim,
            chimeric_MiAIRR     = out_chim,
            chimeric_alignments = out_chim_aln
        )
    end
end

println("CHMMAIRRa fertig.")

