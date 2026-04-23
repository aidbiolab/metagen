import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gmean

# Upload file
df = pd.read_csv('../Raw_Data_for_Figures/Data_for_Figure_2_and_3_merged_rel_abund.csv')
# Remove taxonomy columns and human reads
df = df.drop(columns=['taxonomy_id', 'taxonomy_lvl'], errors='ignore')
df = df[~df['name'].str.contains('Homo sapiens', na=False)]

# ------------------- Define groups -------------------
def get_group(sample):
    if sample.endswith('CP'): return 'Case Plaque (CP)'
    if sample.endswith('CS'): return 'Case Saliva (CS)'
    if sample.endswith('MP'): return 'Patient Plaque (MP)'
    if sample.endswith('MS'): return 'Patient Saliva (MS)'
    if sample.endswith('MT'): return 'Patient Tumor (MT)'
    if sample == 'unknown_barcode': return 'Unknown Barcode'
    return 'Other'

sample_to_group = {col: get_group(col) for col in df.columns[1:]}
grouped = df.set_index('name').T.groupby(sample_to_group).mean().T

# ------------------- CLR Normalization -------------------
def clr_transformation(data):
    """Centered Log-Ratio transformation"""
    data = data + 1e-6                    # pseudocount
    gm = gmean(data, axis=1)              # geometric mean per taxon
    clr = np.log(data / gm[:, np.newaxis])
    return clr

clr_data = clr_transformation(grouped.values)
clr_df = pd.DataFrame(clr_data, index=grouped.index, columns=grouped.columns)

# Select Top 30 taxa
N = 30
top_taxa = clr_df.sum(axis=1).nlargest(N).index
heatmap_data = clr_df.loc[top_taxa]

# ------------------- Plot -------------------
sns.set_style("white")
plt.figure(figsize=(15, 12))

sns.heatmap(
    heatmap_data,
    cmap="RdBu_r",           # Diverging palette - excellent for CLR
    linewidths=0.5,
    linecolor='white',
    center=0,                # Important for CLR
    cbar_kws={'label': 'CLR-transformed Mean Abundance'},
    annot=False
)

plt.title('CLR Normalized Abundance by Sample Type', 
          fontsize=16, fontweight='bold', pad=20)

plt.xlabel('Sample Type', fontsize=13)
plt.ylabel('Taxon', fontsize=12)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()

plt.savefig("heatmap_CLR_normalized.png", dpi=600, bbox_inches='tight')
plt.savefig("heatmap_CLR_normalized.pdf", dpi=600, bbox_inches='tight')

print("✅ CLR-normalized heatmap saved!")