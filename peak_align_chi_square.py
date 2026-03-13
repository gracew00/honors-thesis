# filepath: /home/ubuntu/honors_research/scripts/peak_align_chi_square.py
import sys
import pandas as pd
from scipy.stats import chisquare

def run_chi_square(peak_align_path):
    """
    Read a peak-align CDT file, sum counts into +0..+9 bins,
    and run a chi-square goodness-of-fit test vs uniform.

    peak_align_path is a string path to the peak-align CDT file (tab-delimited with "NAME" column).
    """
    # import the peak align output file
    peak_align_df = pd.read_csv(peak_align_path, sep="\t")

    # drop the NAME column
    if "NAME" in peak_align_df.columns:
        peak_align_df_no_name = peak_align_df.drop(columns=["NAME"])
    else:
        peak_align_df_no_name = peak_align_df

    # create empty dictionary for contingency table
    contig_table = {f"+{i}": 0 for i in range(10)}

    # loop through columns in peak_align df -> get sum for each column + add to dictionary
    for col in peak_align_df_no_name.columns:
        # skip non-position columns if any
        try:
            # get last digit of the column name (e.g., "11" -> 1)
            position = int(str(col)[-1])
        except ValueError:
            # column name does not end in a digit; ignore
            continue

        # sum the values in the column and add to the corresponding position
        contig_table[f"+{position}"] += peak_align_df_no_name[col].sum()

    # convert the dictionary to a dataframe
    contig_df = pd.DataFrame([contig_table])

    # expected data (uniform distribution)
    total = contig_df.sum().sum()
    num_categories = len(contig_df.columns)
    expected = [total / num_categories] * num_categories

    # chi-square test
    chi2_stat, p_value = chisquare(contig_df.values[0], expected)

    print(f"File: {peak_align_path}")
    print("Observed counts (per +position):")
    print(contig_df.to_string(index=False))
    print(f"\nTotal: {total}")
    print(f"Expected per category (uniform): {expected[0]}")
    print(f"\nChi-Square Statistic: {chi2_stat}")
    print(f"P-value: {p_value}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python peak_align_chi_square.py <path_to_cdt>")
        sys.exit(1)

    path = sys.argv[1]
    run_chi_square(path)