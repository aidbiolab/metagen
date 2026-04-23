import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from google.colab import files

# Upload your file
print("Upload your beta_pcoa_braycurtis_counts.csv file")
uploaded = files.upload()
filename = list(uploaded.keys())[0]
df = pd.read_csv(filename)

# ------------------- Create Group column based on suffix -------------------
def get_sample_type(sample_name):
    if sample_name.endswith('CP'):
        return 'Case Plaque (CP)'
    elif sample_name.endswith('CS'):
        return 'Case Saliva (CS)'
    elif sample_name.endswith('MP'):
        return 'Patient Plaque (MP)'
    elif sample_name.endswith('MS'):
        return 'Patient Saliva (MS)'
    elif sample_name.endswith('MT'):
        return 'Patient Tumor (MT)'
    elif sample_name == 'unknown_barcode':
        return 'Unknown Barcode'
    else:
        return 'Other'

df['Group'] = df['sample'].apply(get_sample_type)

# ------------------- Professional Colored PCoA Plot -------------------
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 600

plt.figure(figsize=(13, 10))

# Define nice colors for each group
palette = {
    'Case Plaque (CP)': '#E74C3C',      # Red
    'Case Saliva (CS)': '#3498DB',      # Blue
    'Patient Plaque (MP)': '#2ECC71',   # Green
    'Patient Saliva (MS)': '#F1C40F',   # Yellow
    'Patient Tumor (MT)': '#9B59B6',    # Purple
    'Unknown Barcode': '#7F8C8D'        # Gray
}

sns.scatterplot(
    data=df,
    x="PC1",
    y="PC2",
    hue="Group",
    palette=palette,
    s=140,
    alpha=0.85,
    edgecolor="black",
    linewidth=0.8
)

# Add sample labels
for i in range(len(df)):
    plt.text(
        df.PC1.iloc[i] + 0.012, 
        df.PC2.iloc[i] + 0.012, 
        df['sample'].iloc[i], 
        fontsize=8.5,
        ha='left',
        va='bottom',
        alpha=0.75
    )

plt.title('PCoA Ordination - Bray-Curtis Distance', 
          fontsize=16, fontweight='bold', pad=20)

plt.xlabel('PC1', fontsize=13)
plt.ylabel('PC2', fontsize=13)

# Improve legend
plt.legend(title='Sample Type', title_fontsize=12, fontsize=10, 
           bbox_to_anchor=(1.05, 1), loc='upper left')

plt.grid(True, linestyle='--', alpha=0.5)

# Save
plt.savefig("pcoa_braycurtis_colored_by_type.png", dpi=600, bbox_inches='tight')
plt.savefig("pcoa_braycurtis_colored_by_type.pdf", dpi=600, bbox_inches='tight')

print("✅ Colored PCoA plot saved!")
files.download("pcoa_braycurtis_colored_by_type.png")
files.download("pcoa_braycurtis_colored_by_type.pdf")