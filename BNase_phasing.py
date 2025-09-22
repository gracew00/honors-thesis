import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

def BNase_phasing(tag_path, nucleosome_path):
    tagpileup_output = pd.read_csv(tag_path, sep='\t')
    tagpileup_output = tagpileup_output.drop(columns=['YORF'])

    nucleosome = pd.read_csv(nucleosome_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

    # Getting rid of NAME column in tagpileup
    tagpileup_count = tagpileup_output.drop(columns=['NAME'])

    # Accessing the column "0" in tagpileup_output
    column_0 = tagpileup_output["0"]
    column_0

    # Define the total number of phases (e.g., 10 bins for 0-9, 10-19, ..., 240-249)
    num_phases = 10

    # Initialize a list to store the phase counts for each nucleosome
    nucleosome_phase_counts = []

    # Iterate through each nucleosome (row) in tagpileup_count
    for _, row in tagpileup_output.iterrows():
        # Initialize a dictionary to store phase counts for the current nucleosome
        phase_counts = {}
        
        # Loop through each phase
        for phase in range(num_phases):
            # Get the column indices for the current phase (e.g., 0, 10, 20, ..., for phase 0)
            bin_columns = [str(i) for i in range(phase, 250, 10) if str(i) in tagpileup_count]  # Convert indices to strings to match column names
            
            # Sum the mutation counts for the selected columns for the current nucleosome
            mutation_count = row[bin_columns].sum()
            
            # Store the mutation count for the current phase
            phase_counts[f'phase_{phase + 1}'] = mutation_count  # Phase numbers start from 1
        
        # Add the nucleosome name and phase counts to the list
        nucleosome_phase_counts.append({
            'NAME': row['NAME'],  # Assuming 'NAME' column identifies the nucleosome
            **phase_counts
        })
    # Convert the results into a DataFrame
    nucleosome_phase_counts_df = pd.DataFrame(nucleosome_phase_counts)
    

    nucleosome_phase_counts_df = nucleosome_phase_counts_df.set_index('NAME')
    nucleosome_phase_counts_df.index.name = ''
    nucleosome_phase_counts_df.columns.name = ''

    

    return nucleosome_phase_counts_df

from sklearn.decomposition import PCA
def pca_analysis(nucleosome_phase_counts_df):
    # Assuming nucleosome_phase_counts_df is already loaded as a DataFrame
    # Step 1: Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(nucleosome_phase_counts_df)

    # Step 2: Perform PCA
    pca = PCA(n_components=3)  # Adjust n_components as needed
    principal_components = pca.fit_transform(scaled_data)

    # Step 3: Create a DataFrame for PCA results
    # Adjust column names to match the number of components
    pca_df = pd.DataFrame(data=principal_components, columns=[f'PC{i+1}' for i in range(principal_components.shape[1])])

    # Step 4: Visualize the PCA results (optional)
    plt.figure(figsize=(8, 6))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.7)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.grid()
    plt.show() #how to do this???

    # Step 5: Explained variance ratio
    print("Explained Variance Ratio:", pca.explained_variance_ratio_)

    return pca_df

from mpl_toolkits.mplot3d import Axes3D
def threed_pca(nucleosome_phase_counts_df):
    # Step 1: Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(nucleosome_phase_counts_df)

    # Step 2: Perform PCA
    pca = PCA(n_components=3)  # Adjust n_components as needed
    principal_components = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(data=principal_components, columns=[f'PC{i+1}' for i in range(principal_components.shape[1])])

  
    # Step 3: Create a 3D scatter plot
    fig = plt.figure(figsize=(15, 18))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    ax.scatter(pca_df['PC1'], pca_df['PC2'], pca_df['PC3'], alpha=0.7)

    # Add labels
    ax.set_xlabel('PC1')  # Add padding to the label
    ax.set_ylabel('PC2')  # Add padding to the label
    ax.set_zlabel('PC3')   # Add padding to the label

    ax.set_xlim(pca_df['PC1'].min(), pca_df['PC1'].max())
    ax.set_ylim(pca_df['PC2'].min(), pca_df['PC2'].max())
    ax.set_zlim(pca_df['PC3'].min(), pca_df['PC3'].max())

    plt.show()

    

tag_path = '/home/ubuntu/honors_research/PlusOneDyad_SORT-Expression_250bp_BNase-seq_50U-30min_tagpileup_sense.cdt'
nucleosome_path = '/home/ubuntu/honors_research/RefPT-Krebs/PlusOneDyad_SORT-Expression_250bp.bed'