import pandas as pd

def depmap_snp_editing(path):

    mutations = pd.read_csv(path)

    # filtering for snps (SNV or substiutition)
    mutations = mutations[mutations["Variant Type"].isin(["SNV", "substitution"])]

    # taking chromosome + position columns + making a new column (for end position)
    mutations_new = mutations[["Chromosome"]]
    mutations_new["start"] = mutations["Position"] - 1
    mutations_new["end"] = mutations["Position"]
    
    mutations_new.to_csv("/home/ubuntu/honors_research/depmap/K562_mutations.bed", sep="\t", header=False, index=False)




path = "/home/ubuntu/honors_research/depmap/K562_mutations.csv"