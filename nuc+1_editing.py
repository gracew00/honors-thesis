import pandas as pd



def nuc_editing(file_path):
    """
    Editing TSS_all_adj+1Nuc.bed into separate dfs
    """
    df = df = pd.read_csv(file_path, sep="\t")

    df1 = df[["TSS", "Unnamed: 1", "Unnamed: 2", "Unnamed: 3", 'Unnamed: 5', 'TSS informaion, phase to +1 Nuc']]
    df2 = df[["TSS", "NFR", "Unnamed: 9", "Unnamed: 10", 'Unnamed: 12', 'TSS informaion, phase to +1 Nuc']]
    df3 = df[["TSS", "Unnamed: 17", "Unnamed: 18",  "Unnamed: 19",'Unnamed: 21',  'TSS informaion, phase to +1 Nuc']]
    df4 = df[["TSS", "Unnamed: 23", "Unnamed: 24", "Unnamed: 25",'Unnamed: 27', 'TSS informaion, phase to +1 Nuc']]
    

    #df1.to_csv("/home/ubuntu/honors_research/peak_align/edited_TSS_all_adj+1Nuc.bed", sep="\t", index=False, header = False)
    df2.to_csv("/home/ubuntu/honors_research/peak_align/edited_TSS_all_adj+1Nuc_2.bed", sep="\t", index=False, header = False)
    df3.to_csv("/home/ubuntu/honors_research/peak_align/edited_TSS_all_adj+1Nuc_3.bed", sep="\t", index=False, header = False)
    df4.to_csv("/home/ubuntu/honors_research/peak_align/edited_TSS_all_adj+1Nuc_4.bed", sep="\t", index=False, header = False)


file_path = "/home/ubuntu/honors_research/peak_align/TSS_all_adj+1Nuc.bed"
nuc_editing(file_path)
