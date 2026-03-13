import pandas as pd

def ucsc_vcf_to_bed(vcf_file):
    """
    Convert a VCF file from UCSC Database to BED format.

    Parameters:
    vcf_file (str): Path to the input VCF file.

    vcf_file should have column names of [#chrom, chromStart, chromEnd, name, ref, alfCount, alts, etc]

    Returns:
    pd.DataFrame: DataFrame in BED format with columns ['chrom', 'start', 'end', 'ref', 'alt'].
    """
    # Read the VCF file, skipping header lines
    vcf_df = pd.read_csv(vcf_file, on_bad_lines='skip', sep='\t')

    # Editing the vcf_df to match BED format
    bed_df = vcf_df.iloc[:,:7]

    # removing commas in alts column
    bed_df['alts'] = bed_df['alts'].str.replace(',', '', regex=False)

    # filtering out rows where ref or alt has length > 1
    bed_df = bed_df[(bed_df['ref'].str.len() == 1) & (bed_df['alts'].str.len() == 1)]

    # subtracting chromStart and chromEnd by 1
    bed_df['chromStart'] = bed_df['chromStart'] - 1
    bed_df['chromEnd'] = bed_df['chromEnd'] - 1

    # if there are any "_" in chrom, take the first part (up until the first _)
    bed_df['#chrom'] = bed_df['#chrom'].str.split('_').str[0]

    #final df should have chrom, chromstart, chromend, ref, and alts
    bed_df = bed_df[['#chrom', 'chromStart', 'chromEnd', 'ref', 'alts']]  

    return bed_df

bed_df = ucsc_vcf_to_bed("/home/ubuntu/honors_research/ucsc_snps/dbsnp153.tsv")
bed_df.to_csv("/home/ubuntu/honors_research/ucsc_snps/ucsc_snps.bed", index=False, header = False, sep='\t')