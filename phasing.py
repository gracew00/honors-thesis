
# filepath: /home/ubuntu/honors_research/scripts/plus2Nuc_editing.py
import os
import pandas as pd

# Paths

PLUS2_PATH = "/home/ubuntu/honors_research/peak_align/reference_files/raw_reference/+2_Nuc_all.bed"
PLUS1_TSS_PATH = "/home/ubuntu/honors_research/peak_align/reference_files/raw_reference/Adj+1Nuc_TSS_core.bed"

OUT_PLUS2_DIR = "/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References"
OUT_PLUS1_DIR = "/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References"

os.makedirs(OUT_PLUS2_DIR, exist_ok=True)
os.makedirs(OUT_PLUS1_DIR, exist_ok=True)

def main():
    # -- 1. Load data --

    # +2 nucleosomes, no header
    df_plustwo = pd.read_csv(PLUS2_PATH, sep="\t", header=None)

    # +1 nucleosomes with TSS, has header
    df_tss_plusone = pd.read_csv(PLUS1_TSS_PATH, sep="\t")

    # -- 2. PREPARE +1/TSS TABLE (df_tss_plusone_subset) --

    # keep first 16 columns
    df_tss_plusone_subset = df_tss_plusone.iloc[:, 0:16].copy()

    colnames_plusone_tss = [
        "chr_plusone", "start_plusone_adj", "end_plusone_adj", "id_plusone_adj",
        "distance_plusone_adj", "strand_plusone_adj", "score_plusone_adj",
        "context_plusone_adj", "chr_tss", "start_tss", "end_tss", "id_tss",
        "distance_tss", "strand_tss", "score_tss", "context_tss"
    ]
    df_tss_plusone_subset.columns = colnames_plusone_tss

    # sanity check: |start_tss - start_plusone_adj| == distance_plusone_adj
    df_tss_plusone_subset["sanity_check"] = (
        (df_tss_plusone_subset["start_tss"] - df_tss_plusone_subset["start_plusone_adj"]).abs()
        == df_tss_plusone_subset["distance_plusone_adj"]
    )
    failed = df_tss_plusone_subset[~df_tss_plusone_subset["sanity_check"]]
    print(f"Sanity check failures for +1/TSS: {len(failed)}")

    # -- 3. PREPARE +2 TABLE (df_plustwo) --

    colnames_plustwo = [
        "chr_plustwo", "start_plustwo", "end_plustwo", "id_plustwo",
        "distance_plustwo", "strand_plustwo", "score_plustwo", "context_plustwo",
        "chr_plusone", "start_plusone", "end_plusone", "id_plusone",
        "score_plusone", "strand_plusone"
    ]
    df_plustwo.columns = colnames_plustwo

    # -- 4. MATCH EACH +2 NUC TO NEAREST +1 (BY start_plusone) --

    updated_rows = []
    plusone_starts = df_tss_plusone_subset["start_plusone_adj"]

    for _, row in df_plustwo.iterrows():
        # index of closest +1 by start coordinate
        idx_closest = (plusone_starts - row["start_plusone"]).abs().idxmin()
        closest_row = df_tss_plusone_subset.loc[idx_closest]

        combined_row = row.to_dict()
        for col in df_tss_plusone_subset.columns:
            if col not in df_plustwo.columns:
                combined_row[col] = closest_row[col]

        updated_rows.append(combined_row)

    df_plustwo_updated = pd.DataFrame(updated_rows)

    # -- 5. PHASING FOR +2 NUCLEOSOMES --

    # distance between +2 dyad and TSS
    df_plustwo_updated["distance_tss_plustwo"] = (
        (df_plustwo_updated["start_tss"] - df_plustwo_updated["start_plustwo"]).abs()
    )

    # phasing column
    df_plustwo_updated["phasing_plustwo"] = df_plustwo_updated["distance_tss_plustwo"] % 10

    # drop +1-related columns so only +2 and TSS remain (as in notebook)
    df_plustwo_finalized = df_plustwo_updated.drop(
        columns=[
            "chr_plusone", "start_plusone", "end_plusone",
            "id_plusone", "score_plusone", "strand_plusone",
            "start_plusone_adj", "end_plusone_adj", "id_plusone_adj",
            "distance_plusone_adj", "strand_plusone_adj",
            "score_plusone_adj", "context_plusone_adj"
        ]
    )

    # -- 6. SPLIT +2 BY phasing_plustwo AND WRITE BED FILES --

    for phase in range(10):
        phase_df = df_plustwo_finalized[df_plustwo_finalized["phasing_plustwo"] == phase]
        out_path = os.path.join(OUT_PLUS2_DIR, f"+2_phase{phase}.bed")
        phase_df.to_csv(out_path, sep="\t", index=False, header=False)
        print(f"Wrote +2 phase {phase}: {len(phase_df)} rows -> {out_path}")

    # --7. PHASING FOR +1 NUCLEOSOMES AND WRITE BED FILES --

    # work on a copy to avoid SettingWithCopy warnings
    df_plusone = df_tss_plusone_subset.copy()

    df_plusone["phasing_plusone"] = df_plusone["distance_plusone_adj"] % 10

    # drop TSS-related columns + sanity_check
    df_plusone = df_plusone.drop(
        columns=[
            "chr_tss", "start_tss", "end_tss", "id_tss",
            "distance_tss", "strand_tss", "score_tss",
            "context_tss", "sanity_check"
        ]
    )

    for phase in range(10):
        phase_df = df_plusone[df_plusone["phasing_plusone"] == phase]
        out_path = os.path.join(OUT_PLUS1_DIR, f"+1_phase{phase}.bed")
        phase_df.to_csv(out_path, sep="\t", index=False, header=False)
        print(f"Wrote +1 phase {phase}: {len(phase_df)} rows -> {out_path}")

    print("Done.")

if __name__ == "__main__":
    main()