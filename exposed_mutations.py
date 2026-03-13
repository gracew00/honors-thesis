import pandas as pd
import pysam
import os

def merge_phases(nuc_path, phase0_path, phase1_path, phase2_path, phase3_path, phase4_path, phase5_path, phase6_path, phase7_path, phase8_path, phase9_path, phase10_path):
    # Read each phase file into a DataFrame
    phases = []
    for path in [phase0_path, phase1_path, phase2_path, phase3_path, phase4_path, 
                 phase5_path, phase6_path, phase7_path, phase8_path, 
                 phase9_path, phase10_path]:
        df = pd.read_csv(path, sep="\t", header=None)
        phases.append(df)
    
    nucleosomes = pd.read_csv(nuc_path, sep="\t", header=None, names=["chrom", "start", "end", "NAME"])

    # Concatenate all DataFrames into one
    merged_phases = pd.concat(phases)
    merged_phases = pd.merge(merged_phases, nucleosomes[["chrom", "start", "NAME"]], 
                  on="NAME", how="left")
    merged_phases = merged_phases[["chrom", "start", "NAME"] + list(merged_phases.columns[1:-2])]
    merged_phases = merged_phases.drop(columns=["NAME"])
    
    return merged_df


def selecting_columns(merged_df):
    """
    Selecting columns based on criteria: 
    - position of bp occurs every1 10-ish bp (0,1, 9,10, 11, etc)
    - selected columns must have at least one non-zero value in the row (i.e. must have a signal)
    """
    columns_to_select = []
    for col in range(merged_df.shape[1]):
        if (col) % 10 == 0:  # Check if the column index minus 4 is a multiple of 10
            columns_to_select.extend([col+1, col+2, col + 3])
    columns_to_select = [col for col in columns_to_select if 0 <= col < merged_df.shape[1]]
    merged_phases_filtered = merged_df.iloc[:, [0] + columns_to_select]
    merged_phases_filtered = merged_phases_filtered[merged_phases_filtered.iloc[:, 2:].sum(axis=1) != 0] # filter out rows where all selected columns are 0

    return merged_phases_filtered

def classifying_exposed(mut_path, merged_phases_filtered):
    mutations = pd.read_csv(mut_path, sep="\t", header=None, names=["chrom", "start", "end", "mut_info"])
    filtered_mutations = mutations.copy()
    filtered_mutations["exposed"] = "not_exposed"

    accum = 0
    for idx, row in merged_phases_filtered.iterrows():
        start = int(row[1])
        chrom = row[0]

        
        for col in row.index[2:]:
            value = row[col]
            if value != 0:
                added = int(col)
                pos = start + added
                # print(pos)

                # matching
                matching = filtered_mutations[
                    (filtered_mutations[0] == chrom) &
                    (filtered_mutations[1] == pos)
                ]
                
                # if matching -> put "exposed" in a new column in mutations df

                if not matching.empty:
                    filtered_mutations.loc[
                        (filtered_mutations[0] == chrom) & (filtered_mutations[1] == pos), "exposed"
                    ] = "exposed"
                    accum += 1


    #final_df = pd.DataFrame(output_rows)
    print(accum)
    return filtered_mutations

def importing_vcf(folder_path):
    vcf_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.vcf.gz')]
    variant_set = set()
    for vcf_file in vcf_files:
        vcf = pysam.VariantFile(vcf_file)
        for record in vcf:
            variant_old = f"{record.chrom}:{record.pos}:{record.ref}>{','.join(record.alts)}"
            chr_name, pos, mutation = variant_old.split(':')
            ref, alt = mutation.split(">")
            

            pos_start = int(pos) -1

            # Filter out insertions and deletions (ref and alt lengths must be equal)
            if (len(ref) > 1) or  (len(alt) > 1):
                continue

            pos_start = int(pos) - 1
            pos_end = int(pos) + len(ref) - 1
                
                
                

            variant = f"{record.chrom}:{pos_start}:{pos_end}:{record.ref}>{','.join(record.alts)}"

            variant_set.add(variant)


    # create a dataframe (with chrom, pos, ref>alt columns -> 3 in total)
    vcf_df = pd.DataFrame(list(variant_set), columns=["variant"])
    vcf_df["chrom"] = vcf_df["variant"].apply(lambda x: x.split(":")[0])
    vcf_df["start"] = vcf_df["variant"].apply(lambda x: int(x.split(":")[1]))
    vcf_df["end"] = vcf_df["variant"].apply(lambda x: int(x.split(":")[2]))
    vcf_df["ref_alt"] = vcf_df["variant"].apply(lambda x: x.split(":")[3])
    vcf_df = vcf_df[["chrom", "start", "end", "ref_alt"]]
    return vcf_df

def mutation_context(vcf_df, filtered_mutations):
    filtered_mutations.rename(columns={0: "chrom", 1: "start", 2: "end"}, inplace=True)
    merged_df = filtered_mutations.merge(
        vcf_df[["chrom", "start", "end", "ref_alt"]],  # Select only relevant columns from vcf_df
        on=["chrom", "start", "end"],  # Match on chrom, start, and end
        how="left"  # Use a left join to keep all rows from mutations
    )
    merged_df_nodups = merged_df.drop_duplicates(subset=['chrom', 'start', 'end'])


def add_context_column(df, reference_genome_path):
    """
    Adds a 'context' column to the DataFrame with the nucleotide context (before, SNP, after).

    Parameters:
    df (pd.DataFrame): DataFrame containing mutation information with 'chrom', 'start', and 'end' columns (0-based).
    reference_genome_path (str): Path to the reference genome (e.g., hg38.fa.gz).

    Returns:
    pd.DataFrame: Updated DataFrame with a new 'context' column.
    """
    # Load the reference genome
    fasta = pysam.FastaFile(reference_genome_path)

    # Function to fetch context for a single mutation
    def fetch_context(row):
        try:
            chrom = row['chrom']
            pos = row['start'] + 1  # Convert 0-based to 1-based for pysam
            ref = fasta.fetch(chrom, pos - 1, pos)  # Reference nucleotide at the SNP position
            before = fasta.fetch(chrom, pos - 2, pos - 1)  # Nucleotide before
            after = fasta.fetch(chrom, pos, pos + 1)  # Nucleotide after

            # Making everything in caps
            before = before.upper()
            ref = ref.upper()
            after = after.upper()


            return f"{before}{ref}{after}"
        except KeyError:
            # Handle cases where the chromosome is not found in the reference genome
            return "Context not found"

    # Apply the fetch_context function to each row in the DataFrame
    df['context'] = df.apply(fetch_context, axis=1)

    return df

add_context_column(merged_df_nodups, "/home/ubuntu/honors_research/hg38.fa")

# pathways / variables
mut_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/merged_unique_variants.bed"
nuc_path = "/home/ubuntu/honors_research/peak_align/edited_TSS_all_adj+1Nuc_250bp.bed" #expanded bed file

phase0_path= "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase0.cdt"
phase1_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase1.cdt"
phase2_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase2.cdt"
phase3_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase3.cdt"
phase4_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase4.cdt"
phase5_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase5.cdt"
phase6_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase6.cdt"
phase7_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase7.cdt"
phase8_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase8.cdt"
phase9_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase9.cdt"
phase10_path = "/home/ubuntu/honors_research/peak_align/phased_haplotype_muts/phase10.cdt"

vcf_folder_path = '/home/ubuntu/honors_research/vcf_files'


# function running
merged_df = merge_phases(nuc_path, phase0_path, phase1_path, phase2_path, phase3_path, phase4_path,
                        phase5_path, phase6_path, phase7_path, phase8_path,
                        phase9_path, phase10_path)
merged_phases_filtered = selecting_columns(merged_df)
vcf_df = importing_vcf(vcf_folder_path)
filtered_mutations = classifying_exposed(mut_path, merged_phases_filtered)
merged_df_nodups = mutation_context(vcf_df, filtered_mutations)
merged_df = add_context_column(merged_df_nodups, "/home/ubuntu/honors_research/hg38.fa")

# exporting merged_df  -> use for data analysis 
merged_df.to_csv("/home/ubuntu/honors_research/merged_df_final.csv", index=False)