#pip install vcfpy 
import vcfpy
import os

def get_vcf_files_from_folder(folder_path):
    "Returns a list of all .vcf files in the given folder"
    vcf_files = [os.path.join(folder_path, file) for file in os.listdir(folder_path) if file.endswith('.vcf.gz')]
    print("vcf files are" + vcf_files)
    return vcf_files

def load_vcf(filename):
    reader = vcfpy.Reader.from_path(filename)
    variants = set()
    for record in reader:
        # Create a unique identifier for each variant
        var_id = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        variants.add(var_id)
    return variants

def load_gnomad_vcf(filename):
    gnomad_variants = set()
    reader = vcfpy.Reader.from_path(filename)
    for record in reader:
        var_id = (record.CHROM, record.POS, record.REF, tuple(record.ALT))
        gnomad_variants.add(var_id)
    return gnomad_variants


def filter_variants(union_variants, gnomad_variants):
    """Filter variants to find those that are on the same chromosome as gnomAD
    and do not appear in gnomAD."""
    filtered_variants = set()
    g_var = gnomad_variants[0]
    chrom = g_var[0]
    chrom_variants = {var for var in union_variants if var[0] == chrom}
    unique_variants = chrom_variants - gnomad_variants
   
    return unique_variants

def write_vcf(variants, output_filename):
    """Write variants to a VCF file."""
    header = vcfpy.Header()
    writer = vcfpy.Writer.open(output_filename, header)

    for var in variants:
        record = vcfpy.Record(
            CHROM=var[0],
            POS=var[1],
            ID=None,
            REF=var[2],
            ALT=[vcfpy.Substitution(alt) for alt in var[3]],
            QUAL=None,
            FILTER=None,
            INFO={}
        )
        writer.write_record(record)
    
    writer.close()

# def main(vcf_files, gnomad_file):
#     union_variants = set()
    
#     # Load all VCF files
#     for vcf_file in vcf_files:
#         union_variants.update(load_vcf(vcf_file))
    
#     # Load gnomAD variants
#     gnomad_variants = load_gnomad_vcf(gnomad_file)
    
#     # Filter variants
#     filtered_variants = filter_variants(union_variants, gnomad_variants)

#     # Write filtered unique variants to a VCF file
#     output_filename = 'unique_variants.vcf'
#     write_vcf(filtered_variants, output_filename)


# if __name__ == "__main__":
#     folder_path = '/Users/wudawuda/Documents/honors_research/practice/practice_vcfs'
#     vcf_files = get_vcf_files_from_folder(folder_path)  # Replace with your VCF files
#     gnomad_file = 'gnomad.vcf'  # Replace with your gnomAD VCF file
#     main(vcf_files, gnomad_file)
