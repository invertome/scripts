import pandas as pd
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multicomp import MultiComparison
from sklearn.decomposition import PCA
import argparse
import numpy as np
import re

# Function to sanitize transcript names
def sanitize_name(name):
    match = re.match(r'^>(\S+)', name)
    return match.group(1) if match else name

# Extract TPM values from the quant.sf files for the provided transcripts.
def extract_tpm_from_files(transcripts, input_folder, metadata_df):
    sanitized_transcripts = [sanitize_name(t) for t in transcripts]
    dfs = []

    for index, row in metadata_df.iterrows():
        path = os.path.join(input_folder, row['filename'], 'quant.sf')
        df = pd.read_csv(path, sep='\t', usecols=['Name', 'TPM'])
        df['Name'] = df['Name'].apply(sanitize_name)  # Sanitize names in the quant.sf file
        df = df[df['Name'].isin(sanitized_transcripts)]
        df.set_index('Name', inplace=True)
        df.columns = [row['SampleID']]
        dfs.append(df)

    return pd.concat(dfs, axis=1)

# Plotting function for each transcript using seaborn for better aesthetics.
def plot_data(df, transcript, groups_to_compare, tukey_result, output_folder, plot_type='bar'):
    # Grayscale palette
    gray_palette = sns.color_palette("gray", len(groups_to_compare))
    sns.set_palette(gray_palette)

    if plot_type == 'line':
        # Extract numeric time values from stage names for continuous plotting
        df['Time'] = df['Stage'].apply(lambda x: int(x[1:]))
        plt.figure()
        
        # Calculate SEM for each group
        sem = df.groupby('Time')['TPM'].sem().reset_index()
        
        ax = sns.lineplot(x='Time', y='TPM', data=df, marker='o', sort=True, errorbar=None)
        for idx, row in sem.iterrows():
            plt.errorbar(row['Time'], df[df['Time'] == row['Time']]['TPM'].mean(), yerr=row['TPM'], fmt='o', color='black')
        
        plt.title(transcript)
        plt.xlabel('Day')
        plt.ylabel('TPM')

        # Customize x-axis
        ax.set_xticks(sorted(df['Time'].unique()))
        ax.set_xticklabels(sorted(df['Time'].unique()))
    else:
        plt.figure()
        
        # Calculate mean and SEM for each group
        group_means = df.groupby('Stage')['TPM'].mean().reindex(groups_to_compare)
        group_sem = df.groupby('Stage')['TPM'].sem().reindex(groups_to_compare)
        
        ax = sns.barplot(x='Stage', y='TPM', data=df, order=groups_to_compare, errorbar=None)  # Plot in order of groups_to_compare
        ax.errorbar(range(len(group_means)), group_means, yerr=group_sem, fmt='none', c='black', capsize=5)
        
        plt.title(transcript)
        
        # Annotating significance based on Tukey's post-hoc test
        for i, group in enumerate(groups_to_compare):
            if tukey_result.reject[i]:
                ax.text(i, max(df[df['Stage'] == group]['TPM']) + 0.05, '*', ha='center')

    # Save the plot in both PNG and PDF formats
    for fmt in ['png', 'pdf']:
        plt.savefig(os.path.join(output_folder, f"{transcript}_{plot_type}.{fmt}"))
    plt.close()

# Function to perform ANOVA and post-hoc Tukey's HSD test.
def perform_anova(df, transcript, output_folder):
    mod = sm.formula.ols('TPM ~ Stage', data=df).fit()
    aov_table = sm.stats.anova_lm(mod, typ=2)
    mc = MultiComparison(df['TPM'], df['Stage'])
    tukey_result = mc.tukeyhsd()
    # Write ANOVA and Tukey results to a file
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
    # Save the PCA plot in both PNG and PDF formats
    for fmt in ['png', 'pdf']:
        plt.savefig(os.path.join(output_folder, f"PCA_Plot.{fmt}"))
    plt.close()

# Main function that performs the entire workflow.
def main(metadata_file, transcripts_file, groups_file, input_folder, output_folder, plot_type='bar'):
    metadata = pd.read_csv(metadata_file)
    transcripts = pd.read_csv(transcripts_file, header=None).iloc[:,0].tolist()
    groups_to_compare = pd.read_csv(groups_file, header=None).iloc[:,0].tolist()

    metadata = metadata[metadata['Stage'].isin(groups_to_compare)]
    tpm_matrix = extract_tpm_from_files(transcripts, input_folder, metadata)

    sanitized_transcripts = [sanitize_name(t) for t in transcripts]
    
    for transcript in sanitized_transcripts:
        print(f"Processing transcript: {transcript}")
        try:
            df_transcript = tpm_matrix.loc[transcript].reset_index()
            df_transcript.columns = ['SampleID', 'TPM']
            df_transcript = df_transcript.merge(metadata, on='SampleID', how='left')
            
            tukey_result = perform_anova(df_transcript, transcript, output_folder)
            plot_data(df_transcript, transcript, groups_to_compare, tukey_result, output_folder, plot_type)
        except KeyError as e:
            print(f"KeyError: {e} - Transcript not found in TPM matrix.")
            continue
        
    if len(transcripts) > 3:
        perform_pca(tpm_matrix.T, output_folder)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Conduct ANOVA and plot gene expression data.')
    parser.add_argument('-m', '--metadata_file', type=str, required=True, help='Path to the metadata file.')
    parser.add_argument('-t', '--transcripts_file', type=str, required=True, help='Path to the transcripts file.')
    parser.add_argument('-g', '--groups_file', type=str, required=True, help='Path to the groups file.')
    parser.add_argument('-i', '--input_folder', type=str, required=True, help='Path to the input folder.')
    parser.add_argument('-o', '--output_folder', type=str, required=True, help='Path to the output folder.')
    parser.add_argument('-p', '--plot_type', type=str, choices=['bar', 'line'], default='bar', help='Type of plot to create (default: bar).')
    
    args = parser.parse_args()
    main(args.metadata_file, args.transcripts_file, args.groups_file, args.input_folder, args.output_folder, args.plot_type)
