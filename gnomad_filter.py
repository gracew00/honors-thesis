#pip install vcfpy jupyter server list

#look at exomes first (not genomes)
import pysam
import os

def get_vcf_files_from_folder(folder_path):
    "Returns a list of all .vcf files in the given folder"
    vcf_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.vcf.gz')]
    print(vcf_files)
    return vcf_files

def union_variants(vcf_files):
    variant_set = set()
    for vcf_file in vcf_files:
        vcf = pysam.VariantFile(vcf_file)
        for record in vcf:
            variant = f"{record.chrom}:{record.pos}:{record.ref}>{','.join(record.alts)}"
            variant_set.add(variant)
    return variant_set

def load_gnomad_vcf(gnomad_vcf):
    gnomad_set = set()
    vcf_gnomad = pysam.VariantFile(gnomad_vcf)
    for record in vcf_gnomad:
        gnomad_variant = f"{record.chrom}:{record.pos}:{record.ref}>{','.join(record.alts)}"
        gnomad_set.add(gnomad_variant)

    return gnomad_set

def filter_variants(union_variants, gnomad_variants, chrn):
    """Filter variants to find those that are on the same chromosome as gnomAD
    and do not appear in gnomAD."""
    #chrom_variants = {var for var in union_variants if var[0] == chr}
    chrom_variants = []
    unique_variants = []
    for var in union_variants:
        chr_index_start = var.find("chr")
        chr_index_end = var[chr_index_start:].find(":")
        chr_number = var[chr_index_start:chr_index_end]
        if chr_number == chrn:
            chrom_variants.append(var)
    print(str(len(chrom_variants)) + "number of sample variants on chromosome")
    print(type(chrom_variants))
    accum = 0
    for var in chrom_variants:
        if var not in gnomad_variants:
            print(var)
            accum +=1
            unique_variants.append(var)
    print(str(accum) + "accum variable")
    

    return unique_variants

def export_variants_to_bed(variants, output_file):
    """
    Exports a list of variant strings to a BED file.

    Parameters:
    variants (list of str): List of variants, each formatted as 'chr:start:end'.
    output_file (str): Path to the output BED file.

    Example of variant format: 'chr1:11111:1111112'
    """
    with open(output_file, 'w') as bed_file:
        for variant in variants:
            # Split the string into chromosome, start, and end positions
            
            chr_name, pos, mutation = variant.split(':')
            ref, alt = mutation.split(">")
            pos_start = int(pos) -1

            if len(alt)>len(ref): #insertion
                pos_end = int(pos)+len(alt) -1
            else: #substitution or deletion
                pos_end = int(pos) + len(ref) -1
            
            
            # Write to BED file (converting start and end to integers if needed)
            bed_file.write(f"{chr_name}\t{pos_start}\t{pos_end}\n")

    print(f"Variants successfully exported to {output_file}")





def main(vcf_files, gnomad_file):

    union_vars = union_variants(vcf_files)
    print(str(len(union_vars)) + "length of union variants")

    
    # Load gnomAD variants
    gnomad_variants = load_gnomad_vcf(gnomad_file)
    print(str(len(gnomad_variants)) + "length of gnomad")
    
    # Filter variants
    # gnomad_string = str(gnomad_file)
    # chr_index_start = gnomad_string.find("chr")
    # chr_index_end = gnomad_string[chr_index_start:].find(".")
    chr_number = "chr3"
    print(chr_number + "chr number")
    filtered_variants = filter_variants(union_vars, gnomad_variants, chr_number)

    # Write filtered unique variants to a VCF file
    output_filename = 'unique_variants_chr3.bed'
    export_variants_to_bed(filtered_variants, output_filename)
    print("Unique Variants File Done")




folder_path = '/home/ubuntu/honors_research/vcf_files'
vcf_files = get_vcf_files_from_folder(folder_path)
gnomad_file = '/home/ubuntu/honors_research/gnomad.exomes.v4.1.sites.chr3.vcf.bgz'
all_unique_variants = main(vcf_files, gnomad_file)