import pandas as pd
import pysam
import os


# 1. Load Nucleosome Files

#nucleosomes
def load_nuc1(nuc_path):
    """Load a +1 nucleosome BED file with 9 columns, assign column names.
    
    Expected columns: chrom, start, end, name, col5, col6, col7, col8, phase
    nuc_path is a string path to the BED file"""

    colnames = ["chr", "start", "end", "NAME", "col5", "col6", "col7", "col8", "phase"]
    df = pd.read_csv(nuc_path, sep="\t", header=None)
    df.columns = colnames
    return df

def load_nuc2(nuc_path):
    """Load a +2 nucleosome BED file with 17 columns, assign column names.
    
    Expected columns: chrom, start, end, name, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, phase
    nuc_path is a string path to the BED file"""

    colnames = ["chr", "start", "end", "NAME", "col5", "col6", "col7", "col8", "col9", "col10", "col11", "col12", "col13", "col14", "col15", "col16", "col17","phase"]
    df = pd.read_csv(nuc_path, sep="\t", header=None)
    df.columns = colnames
    return df

#plusone nucs
plusone_nuc_phase0 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase0.bed")
plusone_nuc_phase1 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase1.bed")
plusone_nuc_phase2 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase2.bed")
plusone_nuc_phase3 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase3.bed")
plusone_nuc_phase4 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase4.bed")
plusone_nuc_phase5 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase5.bed")
plusone_nuc_phase6 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase6.bed")
plusone_nuc_phase7 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase7.bed")
plusone_nuc_phase8 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase8.bed")
plusone_nuc_phase9 = load_nuc1("/home/ubuntu/honors_research/peak_align/reference_files/+1_Nuc_References/+1_phase9.bed")
plusone_nuc = pd.concat([plusone_nuc_phase0, plusone_nuc_phase1, plusone_nuc_phase2, plusone_nuc_phase3, plusone_nuc_phase4, plusone_nuc_phase5, plusone_nuc_phase6, plusone_nuc_phase7, plusone_nuc_phase8, plusone_nuc_phase9], ignore_index=True)

#plustwo nucs
plustwo_nuc_phase0 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase0.bed")
plustwo_nuc_phase1 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase1.bed")
plustwo_nuc_phase2 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase2.bed")
plustwo_nuc_phase3 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase3.bed")
plustwo_nuc_phase4 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase4.bed")
plustwo_nuc_phase5 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase5.bed")
plustwo_nuc_phase6 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase6.bed")
plustwo_nuc_phase7 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase7.bed")
plustwo_nuc_phase8 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase8.bed")
plustwo_nuc_phase9 = load_nuc2("/home/ubuntu/honors_research/peak_align/reference_files/+2_Nuc_References/+2_phase9.bed")
plustwo_nuc = pd.concat([plustwo_nuc_phase0, plustwo_nuc_phase1, plustwo_nuc_phase2, plustwo_nuc_phase3, plustwo_nuc_phase4, plustwo_nuc_phase5, plustwo_nuc_phase6, plustwo_nuc_phase7, plustwo_nuc_phase8, plustwo_nuc_phase9], ignore_index=True)


# 2. Load reference genome and get 150bp window sequence at the dyad 

#load reference genome
reference_genome = pysam.FastaFile("/home/ubuntu/honors_research/hg38.fa")

# get the 150bp window sequence at the dyad (add this to a new column in the dataframe)
def get_dyad_sequence(df, reference_genome):
    """For each row in the dataframe, calculate the dyad position (start + 125), 
    then fetch the 150bp sequence centered on the dyad from the reference genome.
    
    Assumes df has columns 'chr', 'start', 'end' and reference_genome is a pysam.FastaFile."""
    
    sequences = []
    for index, row in df.iterrows():
        chrom = row["chr"]
        dyad_position = row["start"] + 125
        start = dyad_position - 75
        end = dyad_position + 75
        sequence = reference_genome.fetch(chrom, start, end).upper()
        sequences.append(sequence)
    df["dyad_sequence"] = sequences
    return df

plusone_nuc = get_dyad_sequence(plusone_nuc, reference_genome)
plustwo_nuc = get_dyad_sequence(plustwo_nuc, reference_genome)


#3. Dicucleotide Categories
#getting dinucleotide context
DINUC_CATEGORY_MAP = {
    "WW-WW":           {"WW-WW"},
    "SS-SS":           {"SS-SS"},
    "SW-WS":           {"SW-WS"},
    "WS-SW":           {"WS-SW"},
    "WW-WS + SW-WW":           {"WW-WS", "SW-WW"},
    "WW-SW + WS-WW":           {"WW-SW", "WS-WW"},
    "WW-SS + SS-WW":           {"WW-SS", "SS-WW"},
    "WS-SS + SS-SW":           {"WS-SS", "SS-SW"},
    "SW-SS + SS-WS":           {"SW-SS", "SS-WS"},
    "WS-WS + SW-SW":           {"WS-WS", "SW-SW"},
}

#helper functions to categorize bases
def base_to_WS(b):
    """Map nucleotide to W/S; return None for N/other.
    b is a string"""
    b = b.upper()
    if b in ("A", "T"):
        return "W"
    elif b in ("C", "G"):
        return "S"
    else:
        return None

def dinuc_from_two(b1, b2):
    """Return 'WW','SS','WS','SW' or None if any base is None.
    Both b1 and b2 are strings"""
    if b1 is None or b2 is None:
        return None
    return b1 + b2

def classify_5bp_window(win5):
    """
    Take a 5bp window, look at first two and last two bases,
    return a concrete pattern string like 'WW-SS' or None if ambiguous.

    win5 is a string of length 5
    """
    if len(win5) < 5:
        return None
    w0 = base_to_WS(win5[0])
    w1 = base_to_WS(win5[1])
    w3 = base_to_WS(win5[3])
    w4 = base_to_WS(win5[4])

    left  = dinuc_from_two(w0, w1)
    right = dinuc_from_two(w3, w4)
    if left is None or right is None:
        return None
    return f"{left}-{right}"

def map_pattern_to_category(pattern):
    """Map a concrete pattern (e.g. 'WW-SS') to one of the 8 merged categories.
    """
    for cat, patterns in DINUC_CATEGORY_MAP.items():
        if pattern in patterns:
            return cat
    return None

def classify_nucleosome_sequence(seq150):
    """
    Slide non-overlapping 5bp windows across seq150,
    count matches for each of the 8 merged categories,
    return the category with the highest count.

    seq150 is a string of length 150 (DNA bases)
    """
    if len(seq150) != 150:
        raise ValueError("Expected 150bp sequence")

    # init counts to 0 for each category
    counts = {cat: 0 for cat in DINUC_CATEGORY_MAP.keys()}

    for i in range(0, 150, 5):
        win = seq150[i:i+5]
        if len(win) < 5:
            continue
        pattern = classify_5bp_window(win)
        if pattern is None:
            continue
        cat = map_pattern_to_category(pattern)
        if cat is not None:
            counts[cat] += 1

    # if everything stayed 0:
    if all(v == 0 for v in counts.values()):
        return "UNKNOWN"

    # pick key with max value
    best_cat = max(counts, key=lambda k: counts[k])
    return best_cat


plusone_nuc["dinuc_category"] = plusone_nuc["dyad_sequence"].apply(classify_nucleosome_sequence)
plustwo_nuc["dinuc_category"] = plustwo_nuc["dyad_sequence"].apply(classify_nucleosome_sequence)


# 4. Exporting
# exporting
def sanitize_category_name(cat):
    """
    Turn category name into something safe for directory/file names.
    """
    return (
        cat.replace(" ", "_")
           .replace("+", "plus")
           .replace("/", "_")
           .replace("-", "_")
    )

def write_bed(df, path):
    """
    Write a 3-column BED (chrom, start, end) from a dataframe.
    Assumes columns are named 'chrom', 'start', 'end'.
    """
    df.to_csv(path, sep="\t", header=False, index=False)

base_dir = "/home/ubuntu/honors_research/peak_align/reference_files"

def export_by_dinuc_and_phase(df, root_name: str):
    """
    For a given nucleosome dataframe:
      - create root folder: base_dir/root_name  (e.g. '+1_dinucleotide')
      - for each dinuc_category, create a subfolder
      - inside each category folder, write one BED per phase
    Assumes df has columns: 'chrom', 'start', 'end', 'phase', 'dinuc_category'.
    """
    out_root = os.path.join(base_dir, root_name)
    os.makedirs(out_root, exist_ok=True)

    dinuc_categories = df["dinuc_category"].dropna().unique()

    for cat in dinuc_categories:
        cat_df = df[df["dinuc_category"] == cat]
        if cat_df.empty:
            continue

        safe_cat = sanitize_category_name(cat)
        cat_dir = os.path.join(out_root, safe_cat)
        os.makedirs(cat_dir, exist_ok=True)

        # loop over phases present in this category
        for phase in sorted(cat_df["phase"].unique()):
            phase_df = cat_df[cat_df["phase"] == phase]
            if phase_df.empty:
                continue

            bed_name = f"{safe_cat}_phase{phase}.bed"
            bed_path = os.path.join(cat_dir, bed_name)
            write_bed(phase_df, bed_path)

# Export +1 and +2 nucleosomes
export_by_dinuc_and_phase(plusone_nuc, "+1_dinucleotide")
export_by_dinuc_and_phase(plustwo_nuc, "+2_dinucleotide")