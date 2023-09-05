# Example command to run the script:
# python select_dge.py [metadata_file] [transcripts_file] [groups_file] [input_folder] [output_folder]

import pandas as pd
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multicomp import MultiComparison
from sklearn.decomposition import PCA

# Extract TPM values from the quant.sf files for the provided transcripts.
def extract_tpm_from_files(transcripts, input_folder, metadata_df):
    dfs = []

    for index, row in metadata_df.iterrows():
        path = os.path.join(input_folder, row['filename'], 'quant.sf')
        df = pd.read_csv(path, sep='\t', usecols=['Name', 'TPM'])
        df = df[df['Name'].isin(transcripts)]
        df.set_index('Name', inplace=True)
        df.columns = [row['SampleID']]
        dfs.append(df)

    return pd.concat(dfs, axis=1)


# Plotting function for each transcript using seaborn for better aesthetics.
def plot_data(df, transcript, groups_to_compare, tukey_result, output_folder):
    # Grayscale palette
    gray_palette = sns.color_palette("gray", len(groups_to_compare))
    sns.set_palette(gray_palette)

    for plot_type in ['bar']:  # Commented out 'box' and 'violin'
        plt.figure()
        ax = sns.barplot(x='Stage', y='TPM', data=df, order=groups_to_compare) # Plot in order of groups_to_compare
        plt.title(transcript)
        
        # Annotating significance based on Tukey's post-hoc test
        for i, group in enumerate(groups_to_compare):
            if tukey_result.reject[i]:
                ax.text(i, max(df[df['Stage'] == group]['TPM']) + 0.05, '*', ha='center')

        for fmt in ['png', 'pdf']:
            plt.savefig(os.path.join(output_folder, f"{transcript}_{plot_type}.{fmt}"))
        plt.close()
        
# Function to perform ANOVA and post-hoc Tukey's HSD test.
def perform_anova(df, transcript, output_folder):
    mod = sm.formula.ols('TPM ~ Stage', data=df).fit()
    aov_table = sm.stats.anova_lm(mod, typ=2)
    mc = MultiComparison(df['TPM'], df['Stage'])
    tukey_result = mc.tukeyhsd()
    with open(os.path.join(output_folder, f"{transcript}_pvalues.txt"), 'w') as f:
        f.write(str(aov_table))
        f.write('\n\n')
        f.write(str(tukey_result))
    return tukey_result  # Return this for later use in plotting

# Function to perform PCA when more than 3 transcripts are provided.
def perform_pca(df, output_folder):
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(df)
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])
    principalDf['Stage'] = df.index
    plt.figure()
    sns.scatterplot(x='principal component 1', y='principal component 2', hue='Stage', data=principalDf)
    for fmt in ['png', 'pdf']:
        plt.savefig(os.path.join(output_folder, f"PCA_Plot.{fmt}"))
    plt.close()

# Main function that performs the entire workflow.
def main(metadata_file, transcripts_file, groups_file, input_folder, output_folder):
    metadata = pd.read_csv(metadata_file)
    transcripts = pd.read_csv(transcripts_file, header=None).iloc[:,0].tolist()
    groups_to_compare = pd.read_csv(groups_file, header=None).iloc[:,0].tolist()

    metadata = metadata[metadata['Stage'].isin(groups_to_compare)]
    tpm_matrix = extract_tpm_from_files(transcripts, input_folder, metadata)

    for transcript in transcripts:
        df_transcript = tpm_matrix.loc[transcript].reset_index()
        df_transcript.columns = ['SampleID', 'TPM']
        df_transcript = df_transcript.merge(metadata, on='SampleID', how='left')
        
        tukey_result = perform_anova(df_transcript, transcript, output_folder)
        plot_data(df_transcript, transcript, groups_to_compare, tukey_result, output_folder)
        
    if len(transcripts) > 3:
        perform_pca(tpm_matrix.T, output_folder)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
