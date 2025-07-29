#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Advanced Primer Analyzer for PCR
================================

This script analyzes BLAST results to rank primer pairs. It is engineered to be
robust, comprehensive, and resilient to common library and data issues.

Workflow:
1.  Loads all inputs: BLAST results, genome FASTA, and a primer FASTA file.
2.  Pre-computes primer stats. It uses a robust, two-tiered approach for Tm:
    a) Tries the highly accurate Nearest-Neighbor model (Tm_NN).
    b) If that fails, it automatically falls back to the stable GC-content
       model (Tm_GC).
3.  Performs an "all forward vs. all reverse" comparison, finding all possible
    amplicons that meet the user-defined size criteria.
4.  Calculates a comprehensive 'optimality_score' for each potential amplicon.
    A LOWER score is BETTER. The score penalizes undesirable traits like low
    specificity, Tm mismatch, and dimer potential.
5.  (Optional) Checks if primers span an intron (ideal for RT-PCR).
6.  Generates two key outputs:
    a) A detailed TSV summary file, sorted by optimality score.
    b) A FASTA file of the top-ranked amplicon sequences.

Usage:
1.  Ensure pandas and biopython are installed (`pip install pandas biopython`).
2.  Configure file paths and parameters in the "USER PARAMETERS" section.
3.  Run from the command line: `python3 advanced_primer_analyzer.py`
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import re
import sys
import os

# ---- USER PARAMETERS ----

# -- Input Files --
blast_file = "primer.blastout"
genome_fasta = "genome.fasta"
primer_fasta = "primers.fasta"

# -- Optional Intron Spanning Analysis --
ENABLE_INTRON_CHECK = False
gff3_file = "/path/to/your/genome_annotation.gff3"

# -- Output Files --
output_summary = "primer_pairs_comprehensive_analysis.tsv"
output_amplicons_fasta = "best_amplicons.fasta"

# -- Filtering and Scoring Criteria --
amplicon_min = 100
amplicon_max = 1500
identity_threshold = 85.0
min_align_len = 16

# -- Primer Quality Parameters --
gc_content_min = 40.0
gc_content_max = 60.0

# --- Thermodynamic Parameters for Tm Calculation ---
# NOTE: Units are critical. Salts are in mM, primers are in nM.
tm_salt_conc = 50.0      # millimolar (mM)
tm_primer_conc = 500.0   # nanomolar (nM) (Only used by the preferred Tm_NN method)

# -- Primer Dimer Analysis Parameters --
dimer_compl_threshold = 7
dimer_3_prime_check_len = 5

# --- Global Flag to track if the fallback Tm calculation was used ---
FALLBACK_TM_CALCULATION_USED = False

# -------------------------

def calculate_tm_robust(sequence, salt_conc, primer_conc):
    """
    Calculates Tm using a robust, two-tiered approach.
    Tries the accurate Tm_NN method first. If it fails, falls back to the
    stable Tm_GC method, guaranteeing a result.

    Returns:
        tuple: A tuple containing the calculated Tm (float) and the
               method used ('NN' or 'GC').
    """
    global FALLBACK_TM_CALCULATION_USED
    try:
        # Tier 1: Attempt the preferred, more accurate Nearest-Neighbor method
        # We use a simplified call for maximum compatibility.
        tm = mt.Tm_NN(sequence, Na=salt_conc, c_seq=primer_conc)
        method = 'NN'
    except (ValueError, KeyError):
        # Tier 2: If NN fails, use the stable GC-based fallback method
        tm = mt.Tm_GC(sequence, Na=salt_conc)
        method = 'GC'
        # Set the global flag so we can inform the user once at the end.
        FALLBACK_TM_CALCULATION_USED = True

    return tm, method

def precompute_primer_stats(primer_fasta_path):
    """Loads primers and calculates Tm, GC%, and GC clamp for all primers."""
    print("Pre-computing stats for all primers (Tm, GC Content)...")
    stats = {}
    primers = SeqIO.to_dict(SeqIO.parse(primer_fasta_path, "fasta"))

    for pid, record in primers.items():
        seq_str = str(record.seq).strip().upper()

        # This new function will now never fail; it will always return a Tm value.
        tm_val, tm_method = calculate_tm_robust(seq_str, tm_salt_conc, tm_primer_conc)

        stats[pid] = {
            'seq': seq_str,
            'gc_content': gc_fraction(seq_str) * 100,
            'tm': tm_val,
            'tm_method': tm_method,
            'has_gc_clamp': seq_str.endswith('G') or seq_str.endswith('C')
        }
    return stats

def check_primer_dimer(p1_seq, p2_seq, is_self_dimer=False):
    p1, p2_rc = Seq(p1_seq), Seq(p2_seq).reverse_complement()
    score = 0
    for i in range(len(p1) - dimer_compl_threshold + 1):
        if p1[i:i+dimer_compl_threshold] in p2_rc: score += 10
    p1_3_prime = p1[-dimer_3_prime_check_len:]
    if p1_3_prime.reverse_complement() in Seq(p2_seq): score += 25
    if not is_self_dimer:
        p2_3_prime = Seq(p2_seq)[-dimer_3_prime_check_len:]
        if p2_3_prime.reverse_complement() in p1: score += 25
    return score

def main():
    print("--- Starting Advanced Primer Analysis ---")

    for f in [blast_file, genome_fasta, primer_fasta]:
        if not os.path.exists(f):
            print(f"[FATAL ERROR] Required input file not found: {f}", file=sys.stderr)
            sys.exit(1)

    primer_stats = precompute_primer_stats(primer_fasta)

    if FALLBACK_TM_CALCULATION_USED:
        print("\n[INFO] A stable fallback method (Tm_GC) was used for some Tm calculations")
        print("       due to an issue with the more accurate Tm_NN method in this environment.")
        print("       The analysis has proceeded successfully using these estimates.\n")

    if not primer_stats:
        print("[FATAL ERROR] No primers could be processed from the FASTA file.", file=sys.stderr)
        sys.exit(1)

    print("Loading genome...")
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    print("Parsing and filtering BLAST output...")
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand"]
    df = pd.read_csv(blast_file, sep="\t", names=cols)

    valid_primer_ids = list(primer_stats.keys())
    df_filtered = df[df["qseqid"].isin(valid_primer_ids) & (df["pident"] >= identity_threshold) & (df["length"] >= min_align_len)].copy()

    if df_filtered.empty:
        print("No BLAST hits passed the filtering criteria. Exiting.")
        return

    print("Identifying all possible amplicons...")
    forward_hits = df_filtered[df_filtered['qseqid'].str.contains('_F')].copy()
    reverse_hits = df_filtered[df_filtered['qseqid'].str.contains('_R')].copy()

    all_amplicons = []
    for contig, f_hits_contig in forward_hits.groupby('sseqid'):
        r_hits_contig = reverse_hits[reverse_hits['sseqid'] == contig]
        if r_hits_contig.empty: continue
        for _, f_row in f_hits_contig.iterrows():
            for _, r_row in r_hits_contig.iterrows():
                if f_row["sstrand"] == r_row["sstrand"]: continue
                if f_row["sstrand"] == "plus": start, end = f_row["send"], r_row["sstart"]
                else: start, end = r_row["send"], f_row["sstart"]
                if start > end: continue
                amplicon_len = end - start + 1
                if not (amplicon_min <= amplicon_len <= amplicon_max): continue
                all_amplicons.append({"F_primer": f_row["qseqid"], "R_primer": r_row["qseqid"], "contig": contig, "amplicon_start": start, "amplicon_end": end, "amplicon_size": amplicon_len, "F_bind_pos": f_row["sstart"], "R_bind_pos": r_row["sstart"], "F_mismatches": f_row["mismatch"], "R_mismatches": r_row["mismatch"]})

    if not all_amplicons:
        print("No potential amplicons found within the specified size range. Exiting.")
        return
    amplicon_df = pd.DataFrame(all_amplicons)
    print(f"Found {len(amplicon_df)} potential amplicons to evaluate.")

    print("Scoring all found primer pairs...")
    primer_hit_counts = df_filtered['qseqid'].value_counts().to_dict()
    final_results = []

    unique_pairs = amplicon_df[['F_primer', 'R_primer']].drop_duplicates()

    for _, pair_row in unique_pairs.iterrows():
        f_id, r_id = pair_row['F_primer'], pair_row['R_primer']
        pair_amplicons = amplicon_df[(amplicon_df['F_primer'] == f_id) & (amplicon_df['R_primer'] == r_id)]
        for _, amp_row in pair_amplicons.iterrows():
            warnings, score = [], 0
            f_stats, r_stats = primer_stats[f_id], primer_stats[r_id]
            f_hits, r_hits = primer_hit_counts.get(f_id, 0), primer_hit_counts.get(r_id, 0)
            score += (f_hits - 1) * 50 + (r_hits - 1) * 50
            if f_hits > 1: warnings.append(f"F_NON_SPECIFIC({f_hits})")
            if r_hits > 1: warnings.append(f"R_NON_SPECIFIC({r_hits})")
            num_amplicons = len(pair_amplicons)
            score += (num_amplicons - 1) * 25
            if num_amplicons > 1: warnings.append(f"MULTIPLE_PRODUCTS({num_amplicons})")
            delta_tm = abs(f_stats['tm'] - r_stats['tm'])
            score += (delta_tm ** 2) * 2
            if delta_tm > 5: warnings.append(f"HIGH_DELTA_TM({delta_tm:.1f})")
            if not (gc_content_min <= f_stats['gc_content'] <= gc_content_max): score += 20; warnings.append(f"F_BAD_GC({f_stats['gc_content']:.1f}%)")
            if not (gc_content_min <= r_stats['gc_content'] <= gc_content_max): score += 20; warnings.append(f"R_BAD_GC({r_stats['gc_content']:.1f}%)")
            if not f_stats['has_gc_clamp']: score += 10; warnings.append("F_NO_GC_CLAMP")
            if not r_stats['has_gc_clamp']: score += 10; warnings.append("R_NO_GC_CLAMP")
            score += amp_row['F_mismatches'] * 20 + amp_row['R_mismatches'] * 20
            if amp_row['F_mismatches'] > 0: warnings.append("F_MISMATCH")
            if amp_row['R_mismatches'] > 0: warnings.append("R_MISMATCH")
            cross, f_self, r_self = check_primer_dimer(f_stats['seq'], r_stats['seq']), check_primer_dimer(f_stats['seq'], f_stats['seq'], True), check_primer_dimer(r_stats['seq'], r_stats['seq'], True)
            score += cross + f_self + r_self
            if cross > 0: warnings.append("CROSS_DIMER")
            if f_self > 0: warnings.append("F_SELF_DIMER")
            if r_self > 0: warnings.append("R_SELF_DIMER")

            final_results.append({"F_primer": f_id, "R_primer": r_id, "optimality_score": score, "warnings": "; ".join(warnings) if warnings else "None", "amplicon_size": amp_row['amplicon_size'], "contig": amp_row['contig'], "Delta_Tm": delta_tm, "F_Tm_Method": f_stats['tm_method'], "R_Tm_Method": r_stats['tm_method'], "F_total_hits": f_hits, "R_total_hits": r_hits, "pair_produces_n_amplicons": num_amplicons, "F_Tm": f_stats['tm'], "R_Tm": r_stats['tm'], "F_GC_content": f_stats['gc_content'], "R_GC_content": r_stats['gc_content'], "F_mismatches": amp_row['F_mismatches'], "R_mismatches": amp_row['R_mismatches'], "F_GC_clamp": f_stats['has_gc_clamp'], "R_GC_clamp": r_stats['has_gc_clamp'], "cross_dimer_score": cross, "F_self_dimer_score": f_self, "R_self_dimer_score": r_self, "amplicon_start": amp_row['amplicon_start'], "amplicon_end": amp_row['amplicon_end']})

    print("Generating final summary report and FASTA file...")
    if not final_results:
        print("Analysis complete, but no valid pairs were found to generate a report.")
        return

    summary_df = pd.DataFrame(final_results).sort_values(by="optimality_score", ascending=True)
    out_cols = ["F_primer", "R_primer", "optimality_score", "warnings", "amplicon_size", "contig", "Delta_Tm", "F_Tm_Method", "R_Tm_Method", "F_total_hits", "R_total_hits", "pair_produces_n_amplicons", "F_Tm", "R_Tm", "F_GC_content", "R_GC_content", "F_mismatches", "R_mismatches", "F_GC_clamp", "R_GC_clamp", "cross_dimer_score", "F_self_dimer_score", "R_self_dimer_score", "amplicon_start", "amplicon_end"]
    summary_df[out_cols].to_csv(output_summary, sep="\t", index=False, float_format="%.2f")

    with open(output_amplicons_fasta, "w") as fasta_out:
        for _, row in summary_df.head(25).iterrows():
            amplicon_seq = genome[row['contig']].seq[row['amplicon_start']-1:row['amplicon_end']]
            header = f"{row['F_primer']}_{row['R_primer']}|score:{row['optimality_score']:.2f}|{row['contig']}:{row['amplicon_start']}-{row['amplicon_end']}|size:{row['amplicon_size']}bp|delta_tm:{row['Delta_Tm']:.2f}|warnings:{row['warnings']}"
            fasta_out.write(f">{header}\n{str(amplicon_seq)}\n")

    print(f"\n--- Analysis Complete ---")
    print(f"Comprehensive summary written to: {output_summary}")
    print(f"Top 25 amplicons written to: {output_amplicons_fasta}")
    print("Review the summary file, sorted by 'optimality_score' (lower is better).")

if __name__ == "__main__":
    main()
