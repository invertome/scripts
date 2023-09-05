import pandas as pd
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.multicomp import MultiComparison
from sklearn.decomposition import PCA

def extract_tpm_from_files(transcripts, input_folder, filenames):
    df_list = []
    for file in filenames:
        path = os.path.join(input_folder, file, 'quant.sf')
        df = pd.read_csv(path, sep='\t', usecols=['Name', 'TPM'])
        df = df[df['Name'].isin(transcripts)]
        df = df.set_index('Name')
        df.columns = [file]
        df_list.append(df)
    return pd.concat(df_list, axis=1)

def plot_data(df, transcript, output_folder):
    color_palette = sns.color_palette("colorblind")
    sns.set_palette(color_palette)

    for plot_type in ['box', 'bar', 'violin']:
        plt.figure()
        if plot_type == 'box':
            sns.boxplot(x='Stage', y='TPM', data=df)
        elif plot_type == 'bar':
            sns.barplot(x='Stage', y='TPM', data=df)
        elif plot_type == 'violin':
            sns.violinplot(x='Stage', y='TPM', data=df)
        
        plt.title(transcript)
        for fmt in ['png', 'pdf']:
            plt.savefig(os.path.join(output_folder, f"{transcript}_{plot_type}.{fmt}"))
        plt.close()

def perform_anova(df, transcript, output_folder):
    mod = sm.formula.ols('TPM ~ Stage', data=df).fit()
    aov_table = sm.stats.anova_lm(mod, typ=2)
    mc = MultiComparison(df['TPM'], df['Stage'])
    tukey_result = mc.tukeyhsd()
    with open(os.path.join(output_folder, f"{transcript}_pvalues.txt"), 'w') as f:
        f.write(str(aov_table))
        f.write('\n\n')
        f.write(str(tukey_result))

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

def main(metadata_file, transcripts_file, groups_file, input_folder, output_folder):
    metadata = pd.read_csv(metadata_file)
    transcripts = pd.read_csv(transcripts_file, header=None).iloc[:,0].tolist()
    groups_to_compare = pd.read_csv(groups_file, header=None).iloc[:,0].tolist()

    metadata = metadata[metadata['Stage'].isin(groups_to_compare)]
    tpm_matrix = extract_tpm_from_files(transcripts, input_folder, metadata['filename'].tolist())

    for transcript in transcripts:
        df_transcript = tpm_matrix.loc[transcript].reset_index()
        df_transcript.columns = ['SampleID', 'TPM']
        df_transcript = df_transcript.merge(metadata, on='SampleID', how='left')
        
        plot_data(df_transcript, transcript, output_folder)
        perform_anova(df_transcript, transcript, output_folder)
        
    if len(transcripts) > 3:
        perform_pca(tpm_matrix.T, output_folder)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
