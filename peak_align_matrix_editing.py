import pandas as pd

def peak_align_editing(peak_align_path, reference_bed_path):
    """
    Edits peak alignment matrix to separate phases and prepare for composite plot graphing.

    Parameters:
    peak_align_path (str): Path to the peak alignment matrix file.
    reference_bed_path (str): Path to the reference BED file.
    """

    #Load data
    peak_align = pd.read_csv(peak_align_path, sep="\t")
    reference = pd.read_csv(reference_bed_path, sep="\t", header=None)

    #Extracting necessary columns from reference
    reference.columns = ["chrom", "start", "end", "name", "number", "score", "numb2", "strand"]
    phase = reference["strand"].str.split("_").str[-1]
    reference["phase"] = phase

    #Editing peak alignment matrix
    peak_align = peak_align.drop(columns = "YORF")
    reference_updated = reference[["name", "phase"]]
    
    for i in range(len(reference_updated["phase"])):
        phase_name = reference_updated["phase"][i]
        phase_i = phase_name.find("phase")
        phase = phase_name[phase_i:]
        #print(phase)

        reference_updated.loc[i, "phase"] = phase
    peak_align["phase"] = reference["phase"]

    #separating by phases
    dataframes = {}

    for phase_value, group in peak_align.groupby('phase'):
        dataframes[f"{phase_value}"] = group

    phase0 = dataframes["phase0"].drop(columns = "phase")
    phase1 = dataframes["phase1"].drop(columns = "phase")
    phase2 = dataframes["phase2"].drop(columns = "phase")
    phase3 = dataframes["phase3"].drop(columns = "phase")
    phase4 = dataframes["phase4"].drop(columns = "phase")
    phase5 = dataframes["phase5"].drop(columns = "phase")
    phase6 = dataframes["phase6"].drop(columns = "phase")
    phase7 = dataframes["phase7"].drop(columns = "phase")
    phase8 = dataframes["phase8"].drop(columns = "phase")
    phase9 = dataframes["phase9"].drop(columns = "phase")
    phase10 = dataframes["phaseno10x"].drop(columns = "phase")

    #Checking
    sum = len(phase0) + len(phase1) + len(phase2) + len(phase3) + len(phase4) + len(phase5) + len(phase6) + len(phase7) + len(phase8) + len(phase9) + len(phase10)
    print("Sum: " + str(sum))
    print("Total original length: " + str(len(peak_align)))

    phase0.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase0_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase1.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase1_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase2.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase2_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase3.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase3_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase4.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase4_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase5.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase5_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase6.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase6_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase7.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase7_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase8.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase8_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase9.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase9_ucsc_4_250bp.cdt", sep="\t",index=False)
    phase10.to_csv("/home/ubuntu/honors_research/peak_align/ucsc_muts/phase10_ucsc_4_250bp.cdt", sep="\t",index=False)




peak_align_path = "/home/ubuntu/honors_research/peak_align/ucsc_snps_edited_TSS_all_adj+1Nuc_4_250bp_Output.cdt"
reference_bed = "/home/ubuntu/honors_research/peak_align/edited_TSS_all_adj+1Nuc_4_250bp.bed"
