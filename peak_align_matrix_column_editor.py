import os
import glob
import pandas as pd

def drop_first_column_in_folder(folder, pattern="*.cdt", sep="\t"):
    """
    Access folder of peak aligns and dropping the first column ("YORF")
    from each .cdt file. The modified files will overwrite the original ones.
    """
    folder = os.path.abspath(folder)
    print("Folder:", folder)

    files = glob.glob(os.path.join(folder, pattern))
    print("Found files:", files)

    for path in files:
        print(f"Processing {path}")
        df = pd.read_csv(path, sep=sep, dtype=str)

        if "YORF" not in df.columns:
            print("  No 'YORF' column, skipping")
            continue

        df = df.drop(columns=["YORF"])
        # keep header, no index
        df.to_csv(path, sep=sep, index=False)

drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/SS_SS")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/SW_WS")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WS_SW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WW_SS_plus_SS_WW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WW_WS_plus_SW_WW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/SW_SS_plus_SS_WS")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WS_SS_plus_SS_SW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WS_WS_plus_SW_SW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WW_SW_plus_WS_WW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+2_dinucleotide/WW_WW")

drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/SS_SS")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/SW_WS")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WS_SW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WW_SS_plus_SS_WW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WW_WS_plus_SW_WW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/SW_SS_plus_SS_WS")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WS_SS_plus_SS_SW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WS_WS_plus_SW_SW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WW_SW_plus_WS_WW")
drop_first_column_in_folder("/home/ubuntu/honors_research/peak_align/+1_dinucleotide/WW_WW")