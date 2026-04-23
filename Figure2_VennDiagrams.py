import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
from google.colab import files
import numpy as np

# ------------------- Upload your merged file -------------------
print("Upload your merged_rel_abund.csv file")
uploaded = files.upload()
df = pd.read_csv(list(uploaded.keys())[0])

# Keep only taxon name and samples
df = df[['name'] + [col for col in df.columns if col not in ['taxonomy_id', 'taxonomy_lvl']]]

# ------------------- Define groups -------------------
def get_type(sample):
    if sample.endswith('MP'): return 'Patient Plaque'
    if sample.endswith('MS'): return 'Patient Saliva'
    if sample.endswith('MT'): return 'Patient Tumor'
    if sample.endswith('CP'): return 'Control Plaque'
    if sample.endswith('CS'): return 'Control Saliva'
    return 'Other'

# Create presence/absence (present if relative abundance > 0.001)
abund = df.set_index('name').T
abund = abund > 0.001   # threshold to consider taxon "present"

# ------------------- Create the 4 comparisons -------------------
plt.figure(figsize=(16, 12))

# A: Patient Plaque vs Patient Saliva vs Patient Tumor (3-way Venn)
plt.subplot(2, 2, 1)
pp = set(abund[abund.index.str.endswith('MP')].columns[abund[abund.index.str.endswith('MP')].any()])
ps = set(abund[abund.index.str.endswith('MS')].columns[abund[abund.index.str.endswith('MS')].any()])
pt = set(abund[abund.index.str.endswith('MT')].columns[abund[abund.index.str.endswith('MT')].any()])

venn3([pp, ps, pt], set_labels=('Patient Plaque', 'Patient Saliva', 'Patient Tumor'))
plt.title('A. Patient Plaque vs Patient Saliva vs Patient Tumor', fontsize=13, fontweight='bold')

# B: Control Plaque vs Control Saliva (2-way Venn)
plt.subplot(2, 2, 2)
cp = set(abund[abund.index.str.endswith('CP')].columns[abund[abund.index.str.endswith('CP')].any()])
cs = set(abund[abund.index.str.endswith('CS')].columns[abund[abund.index.str.endswith('CS')].any()])

venn2([cp, cs], set_labels=('Control Plaque', 'Control Saliva'))
plt.title('B. Control Plaque vs Control Saliva', fontsize=13, fontweight='bold')

# C: Control Saliva vs Patient Saliva
plt.subplot(2, 2, 3)
venn2([cs, ps], set_labels=('Control Saliva', 'Patient Saliva'))
plt.title('C. Control Saliva vs Patient Saliva', fontsize=13, fontweight='bold')

# D: Control Plaque vs Patient Plaque
plt.subplot(2, 2, 4)
venn2([cp, pp], set_labels=('Control Plaque', 'Patient Plaque'))
plt.title('D. Control Plaque vs Patient Plaque', fontsize=13, fontweight='bold')

plt.suptitle('Venn Diagrams of Shared and Unique Taxa Among Sample Types\n(Threshold: relative abundance > 0.001)', 
             fontsize=16, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save
plt.savefig("venn_diagrams_4panel.png", dpi=600, bbox_inches='tight')
plt.savefig("venn_diagrams_4panel.pdf", dpi=600, bbox_inches='tight')

print("✅ 4 Venn diagrams saved successfully!")
files.download("venn_diagrams_4panel.png")
files.download("venn_diagrams_4panel.pdf")