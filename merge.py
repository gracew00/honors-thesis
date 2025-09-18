import pandas as pd


def merge_bed(bed_files, output_file):
    """
    Merges multiple BED files into a single BED file.
    """
    dataframes = [pd.read_csv(f, sep="\t", header=None) for f in bed_files]
    merged_df = pd.concat(dataframes)
    merged_df.sort_values(by=[0, 1], inplace=True)
    merged_df.to_csv(output_file, sep="\t", header=False, index=False)
    print(f"Merged BED file saved to {output_file}")






bed_files = ["unique_variants_chr1.bed", "unique_variants_chr2.bed", "unique_variants_chr3.bed", "unique_variants_chr4.bed", 
             "unique_variants_chr5.bed", "unique_variants_chr6.bed", "unique_variants_chr7.bed", "unique_variants_chr8.bed",
             "unique_variants_chr9.bed", "unique_variants_chr10.bed", "unique_variants_chr11.bed", "unique_variants_chr12.bed",
             "unique_variants_chr13.bed", "unique_variants_chr14.bed", "unique_variants_chr15.bed", "unique_variants_chr16.bed",
             "unique_variants_chr17.bed", "unique_variants_chr18.bed", "unique_variants_chr19.bed", "unique_variants_chr20.bed",
             "unique_variants_chr21.bed", "unique_variants_chr22.bed", "unique_variants_chrX.bed", "unique_variants_chrY.bed"]

output_file = "merged_unique_variants.bed"

merge_bed(bed_files, output_file)