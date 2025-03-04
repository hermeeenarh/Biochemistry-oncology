Code 2.4:
Biochemistry & Oncology

Task Code 2.6:
Transcriptomics

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import requests

# Read the transcriptomics dataset from the url, specifying the appropriate delimiter used to separate the columns in the dataset
url = "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"
transcriptomics_data = pd.read_csv(url, sep=r'\s+', engine='python')

# Compute -log10(pvalue) for all genes which forms the vertical axis of the volcano plot
transcriptomics_data['neg_log10_pvalue'] = -np.log10(transcriptomics_data['pvalue'])

# Classify genes into categories (i.e upregulated, downregulated, and not differentially expressed) based on the conditions
transcriptomics_data['category'] = np.where((transcriptomics_data['log2FoldChange'] > 1) & (transcriptomics_data['pvalue'] < 0.01), 'Upregulated',
                                            np.where((transcriptomics_data['log2FoldChange'] < -1) & (transcriptomics_data['pvalue'] < 0.01), 'Downregulated', 'Not significant'))

# Generate the volcano plot with unique colors for each category
plt.figure(figsize=(10, 6))

# Define colors for each category
category_colors = {'Upregulated': 'red', 'Downregulated': 'blue', 'Not significant': 'gray'}

# Create the plot with the three categories
sns.scatterplot(data=transcriptomics_data, x='log2FoldChange', y='neg_log10_pvalue', hue='category', palette=category_colors, edgecolor=None, alpha=0.7)

# Add lines for thresholds (log2FoldChange > 1, log2FoldChange < -1, pvalue < 0.01)
plt.axhline(y=-np.log10(0.01), color='black', linestyle='--')  # Significance threshold for pvalue = 0.01
plt.axvline(x=1, color='green', linestyle='--')  # Threshold for upregulation (log2FoldChange > 1)
plt.axvline(x=-1, color='green', linestyle='--')  # Threshold for downregulation (log2FoldChange < -1)

# Top 5 upregulated genes
top_upregulated_genes = transcriptomics_data[(transcriptomics_data['log2FoldChange'] > 1) & (transcriptomics_data['pvalue'] < 0.01)].nlargest(5, 'log2FoldChange')

# Top 5 downregulated genes
top_downregulated_genes = transcriptomics_data[(transcriptomics_data['log2FoldChange'] < -1) & (transcriptomics_data['pvalue'] < 0.01)].nsmallest(5, 'log2FoldChange')

# Annotate top 5 upregulated genes on the plot (black color for better visibility)
for i, gene in top_upregulated_genes.iterrows():
    plt.text(gene['log2FoldChange'], gene['neg_log10_pvalue'], gene['Gene'], fontsize=9, ha='left', color='black', verticalalignment='bottom', horizontalalignment='right')

# Annotate top 5 downregulated genes on the plot (black color for better visibility)
for i, gene in top_downregulated_genes.iterrows():
    plt.text(gene['log2FoldChange'], gene['neg_log10_pvalue'], gene['Gene'], fontsize=9, ha='right', color='black', verticalalignment='top', horizontalalignment='left')

# Specify plot details icluding the plot title, and the horizontal and vertical axes labels
plt.title('Volcano Plot of Differentially Expressed Genes from RNA-Seq Analysis of X Treatment on Diseased Cell Lines')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 P-value')

# Modify the legend key to show 'Significance'
plt.legend(title='Significance', loc='upper right', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()

# Output the upregulated and downregulated genes
print("Upregulated genes:")
print(top_upregulated_genes[['Gene', 'log2FoldChange', 'pvalue']])

print("\nDownregulated genes:")
print(top_downregulated_genes[['Gene', 'log2FoldChange', 'pvalue']])

# Get the functions of the top 5 upregulated and top 5 downregulated genes from Ensembl (did not use GeneCards as it blocks automated requests from scripts)
def get_gene_function_ensembl(gene_name):
    # Ensembl API URL to get gene information (including GO terms)
    url = f"https://rest.ensembl.org/lookup/symbol/human/{gene_name}?content-type=application/json"
    
    try:
        # Send GET request to Ensembl API
        response = requests.get(url)
        
        if response.status_code == 200:
            gene_info = response.json()
            
            # Get Gene Ontology (GO) terms if available
            go_terms = gene_info.get('go_terms', [])
            if go_terms:
                return ', '.join(go_terms)  # Join multiple GO terms if available
            
            # If GO terms are not available, try the description
            if 'description' in gene_info:
                return gene_info['description']
            else:
                return "No functional information available"
        else:
            return f"Failed to retrieve data (Status code: {response.status_code})"
    except Exception as e:
        return f"Error occurred: {str(e)}"

# Top 5 upregulated genes
print("\nTop 5 Upregulated Genes' Functions:")
for gene in top_upregulated_genes['Gene']:
    function = get_gene_function_ensembl(gene)
    print(f"{gene}: {function}")

# Top 5 downregulated genes
print("\nTop 5 Downregulated Genes' Functions:")
for gene in top_downregulated_genes['Gene']:
    function = get_gene_function_ensembl(gene)
    print(f"{gene}: {function}")

Top 5 Upregulated Genes:
       Gene 
DTHD1
EMILIN2
PI16
C4orf45
FAM180B 

Top 5 Downregulated Genes:
       Gene
TBX5
IFITM1
TNN
COL13A1
IFITM3


