import os
import pandas as pd
import pysam


# 1. Paths/Directories

MUT_BED_PATH = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/merged_unique_variants.bed"

PLUS1_NUC_REF_DIR = "/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References"
PLUS2_NUC_REF_DIR = "/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References"

PLUS1_ALIGN_DIR = "/home/ubuntu/honors_research/peak_align/+1_Nuc_Aligns"
PLUS2_ALIGN_DIR = "/home/ubuntu/honors_research/peak_align/+2_Nuc_Aligns"

VCF_DIR = "/home/ubuntu/honors_research/vcf_files"
FASTA_PATH = "/home/ubuntu/honors_research/hg38.fa"

OUT_PLUS1_CSV = "/home/ubuntu/honors_research/exposed_muts_plusone_trial.csv"
OUT_PLUS2_CSV = "/home/ubuntu/honors_research/exposed_muts_plustwo_trial.csv"


# 2. Helper functions to load data
def load_mutations(mut_path: str) -> pd.DataFrame:
    df = pd.read_csv(mut_path, sep="\t", header=None)
    df.columns = ["chrom", "start", "end"]
    return df


def load_nuc1(nuc_path: str) -> pd.DataFrame:
    colnames = ["chr", "start", "end", "NAME", "col5", "col6", "col7", "col8", "phase"]
    df = pd.read_csv(nuc_path, sep="\t", header=None)
    df.columns = colnames
    return df


def load_nuc2(nuc_path: str) -> pd.DataFrame:
    colnames = [
        "chr", "start", "end", "NAME",
        "col5", "col6", "col7", "col8",
        "col9", "col10", "col11", "col12",
        "col13", "col14", "col15", "col16",
        "col17", "phase",
    ]
    df = pd.read_csv(nuc_path, sep="\t", header=None)
    df.columns = colnames
    return df


def load_all_nuc_refs(base: str, plus: str) -> pd.DataFrame:
    """
    plus: "+1" or "+2"
    """
    dfs = []
    loader = load_nuc1 if plus == "+1" else load_nuc2
    for phase in range(10):
        path = os.path.join(base, f"{plus}_phase{phase}_250bp.bed")
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        dfs.append(loader(path))
    return pd.concat(dfs, ignore_index=True)


def load_all_peak_aligns(base: str, plus: str) -> pd.DataFrame:
    """
    Load all phase CDT peak-align files for +1 or +2.
    """
    dfs = []
    for phase in range(10):
        path = os.path.join(
            base,
            f"gnomAD_{plus}_phase{phase}_250bp_Output.cdt",
        )
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        dfs.append(pd.read_csv(path, sep="\t"))
    return pd.concat(dfs, ignore_index=True)


# 3. Reordering peak-align matrices

def split_name(peak_aligns: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a peak-align CDT matrix into a matrix with:
      - 'chr', 'start' (dyad-125) as first two columns
      - numeric columns (0..249) as before (but NAME, end removed).

    Expects 'NAME' like "chr1_2229073_2229073".
    """
    # split into chr, start, end
    peak_aligns[["chr", "start", "end"]] = peak_aligns["NAME"].str.split("_", expand=True)
    peak_aligns["start"] = peak_aligns["start"].astype(int) - 125
    # drop NAME
    peak_aligns = peak_aligns.drop(columns=["NAME"])
    # reorder: chr, start, then the rest (excluding end)
    reordered_columns = ["chr", "start"] + [
        col for col in peak_aligns.columns if col not in ["chr", "start", "end"]
    ]
    peak_aligns = peak_aligns[reordered_columns]
    return peak_aligns

def load_all_vcfs(vcf_dir: str) -> pd.DataFrame:
    """
    Read all .vcf.gz in vcf_dir and return a DataFrame:
      chrom, start, end, ref_alt
    Using 0-based start / end (bed-like).
    """
    vcf_files = [
        os.path.join(vcf_dir, f)
        for f in os.listdir(vcf_dir)
        if f.endswith(".vcf.gz")
    ]
    variant_set = set()
    for vcf_path in vcf_files:
        vcf = pysam.VariantFile(vcf_path)
        for record in vcf:
            variant_old = f"{record.chrom}:{record.pos}:{record.ref}>{','.join(record.alts)}"
            chr_name, pos_str, mutation = variant_old.split(":")
            ref, alt = mutation.split(">")
            # skip indels
            if len(ref) > 1 or len(alt) > 1:
                continue

            pos = int(pos_str)
            pos_start = pos - 1
            pos_end = pos + len(ref) - 1
            variant = f"{record.chrom}:{pos_start}:{pos_end}:{record.ref}>{','.join(record.alts)}"
            variant_set.add(variant)

    vcf_df = pd.DataFrame(list(variant_set), columns=["variant"])
    vcf_df["chrom"] = vcf_df["variant"].apply(lambda x: x.split(":")[0])
    vcf_df["start"] = vcf_df["variant"].apply(lambda x: int(x.split(":")[1]))
    vcf_df["end"] = vcf_df["variant"].apply(lambda x: int(x.split(":")[2]))
    vcf_df["ref_alt"] = vcf_df["variant"].apply(lambda x: x.split(":")[3])
    vcf_df = vcf_df[["chrom", "start", "end", "ref_alt"]]
    return vcf_df

# 4. Exposed and Unexposed Mutations

def split_exposed_unexposed(merged_phases: pd.DataFrame):
    """
    Using your column-index pattern:
      - build 'columns_to_select' for exposed
      - complement as unexposed
      - filter rows where sum != 0
    Returns (exposed_df, unexposed_df).
    """
    # choose exposed columns based on index pattern
    columns_to_select = []
    ncols = merged_phases.shape[1]
    for col in range(ncols):
        if (col) % 10 == 0:
            columns_to_select.extend([col + 1, col + 2, col + 3])

    columns_to_select = [c for c in columns_to_select if 0 <= c < ncols]

    # exposed: chr (0) + exposed cols
    exposed = merged_phases.iloc[:, [0] + columns_to_select]
    exposed = exposed[exposed.iloc[:, 2:].sum(axis=1) != 0]

    # unexposed: all other numeric columns
    all_columns = set(range(ncols))
    unexposed_columns = list(all_columns - set(columns_to_select))
    unexposed_columns = [c for c in unexposed_columns if 0 <= c < ncols]
    unexposed_columns = unexposed_columns[1:]  # drop column 0 (chr) so we will explicitly add 0,1

    unexposed = merged_phases.iloc[:, [0, 1] + unexposed_columns]
    unexposed = unexposed[unexposed.iloc[:, 2:].sum(axis=1) != 0]

    return exposed, unexposed


# 5. Getting positions from peak align matrix

def positions_from_matrix(df: pd.DataFrame):
    """
    Return set of (chrom, pos) where df has non-zero values.
    Expects first two columns: chr, start, then numeric columns 0..N.
    """
    pos_set = set()
    for _, row in df.iterrows():
        chrom = row.iloc[0]
        start = int(row.iloc[1])
        for col in row.index[2:]:
            value = row[col]
            if value != 0:
                offset = int(col)
                pos = start + offset
                pos_set.add((chrom, pos))
    return pos_set


# 6. Categorizing mutations as exposed/unexposed

def label_mutations_exposure(mutations: pd.DataFrame,
                             exposed_pos: set,
                             unexposed_pos: set) -> pd.DataFrame:
    df = mutations.copy()
    df["exposed"] = "none"
    accum_exposed = 0
    accum_unexposed = 0

    for idx, row in df.iterrows():
        key = (row["chrom"], row["start"])
        if key in exposed_pos:
            df.at[idx, "exposed"] = "exposed"
            accum_exposed += 1
        elif key in unexposed_pos:
            df.at[idx, "exposed"] = "unexposed"
            accum_unexposed += 1

    print("Total exposed:", accum_exposed)
    print("Total unexposed:", accum_unexposed)
    print(df.groupby("exposed").size())
    df = df[df["exposed"] != "none"]
    return df




# 7. Adding trinucleotide context for each mutation

def add_context_column(df: pd.DataFrame, reference_genome_path: str) -> pd.DataFrame:
    """
    Adds 'context' = before + ref + after (3bp) using hg38.
    Assumes df has chrom/start/end (0-based).
    """
    fasta = pysam.FastaFile(reference_genome_path)

    def fetch_context(row):
        try:
            chrom = row["chrom"]
            pos = row["start"] + 1  # 1-based for pysam
            ref = fasta.fetch(chrom, pos - 1, pos)
            before = fasta.fetch(chrom, pos - 2, pos - 1)
            after = fasta.fetch(chrom, pos, pos + 1)
            return f"{before.upper()}{ref.upper()}{after.upper()}"
        except KeyError:
            return "Context not found"

    df = df.copy()
    df["context"] = df.apply(fetch_context, axis=1)
    return df


# 8. Main pipeline

def main():
    # 1) Load mutations
    print("Loading mutations...")
    mutations = load_mutations(MUT_BED_PATH)

    # 2) Load peak-align matrices
    print("Loading +1 peak-align matrices...")
    plus1_peak = load_all_peak_aligns(PLUS1_ALIGN_DIR, "+1")
    print("Loading +2 peak-align matrices...")
    plus2_peak = load_all_peak_aligns(PLUS2_ALIGN_DIR, "+2")

    # 3) Convert matrices to chr/start anchored
    print("Computing dyad-anchored positions...")
    plus1_phases = split_name(plus1_peak)
    plus2_phases = split_name(plus2_peak)

    # 4) Split into exposed/unexposed matrices
    print("Splitting +1 into exposed/unexposed...")
    exposed_plus1, unexposed_plus1 = split_exposed_unexposed(plus1_phases)
    print("Splitting +2 into exposed/unexposed...")
    exposed_plus2, unexposed_plus2 = split_exposed_unexposed(plus2_phases)

    # 5) Build position sets
    print("Building position sets...")
    exposed_pos_plus1 = positions_from_matrix(exposed_plus1)
    unexposed_pos_plus1 = positions_from_matrix(unexposed_plus1)
    exposed_pos_plus2 = positions_from_matrix(exposed_plus2)
    unexposed_pos_plus2 = positions_from_matrix(unexposed_plus2)

    # 6) Label mutations for +1 and +2
    print("Labeling mutations for +1...")
    plus1_filtered_mutations = label_mutations_exposure(
        mutations, exposed_pos_plus1, unexposed_pos_plus1
    )
    print("Labeling mutations for +2...")
    plus2_filtered_mutations = label_mutations_exposure(
        mutations, exposed_pos_plus2, unexposed_pos_plus2
    )

    # 7) Load VCF variants and merge
    print("Loading VCF variants...")
    vcf_df = load_all_vcfs(VCF_DIR)

    print("Merging +1 exposure with variants...")
    merged_plus1 = plus1_filtered_mutations.merge(
        vcf_df[["chrom", "start", "end", "ref_alt"]],
        on=["chrom", "start", "end"],
        how="left",
    )

    print("Merging +2 exposure with variants...")
    merged_plus2 = plus2_filtered_mutations.merge(
        vcf_df[["chrom", "start", "end", "ref_alt"]],
        on=["chrom", "start", "end"],
        how="left",
    )

    # remove duplicates by coordinate
    merged_plus1 = merged_plus1.drop_duplicates(subset=["chrom", "start", "end"])
    merged_plus2 = merged_plus2.drop_duplicates(subset=["chrom", "start", "end"])

    # 8) Add 3bp context
    print("Adding 3bp context for +1...")
    merged_plus1 = add_context_column(merged_plus1, FASTA_PATH)
    print("Adding 3bp context for +2...")
    merged_plus2 = add_context_column(merged_plus2, FASTA_PATH)

    # 9) Export
    print(f"Writing {OUT_PLUS1_CSV}")
    merged_plus1.to_csv(OUT_PLUS1_CSV, index=False)
    print(f"Writing {OUT_PLUS2_CSV}")
    merged_plus2.to_csv(OUT_PLUS2_CSV, index=False)
    print("Done.")


if __name__ == "__main__":
    main()