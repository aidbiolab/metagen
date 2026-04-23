
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your data
#df = pd.read_csv("alpha_diversity_counts.csv")
df = pd.read_csv('../Raw_Data_for_Figures/Data_for_Figure_4a_alpha_diversity_counts.csv')

# Use default fonts available in Colab + clean style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 600

# Create the figure
fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)
fig.suptitle('Alpha Diversity Metrics', 
             fontsize=18, fontweight='bold', y=1.02)

metrics = ["shannon", "observed_features", "chao1", "berger_parker_d"]
titles = ["Shannon Diversity", "Observed Features (Richness)", 
          "Chao1 Richness Estimator", "Berger-Parker Dominance"]
ylabels = ["Shannon Index", "Number of Species", 
           "Estimated Richness", "Dominance Index"]

for i, (metric, title, ylabel) in enumerate(zip(metrics, titles, ylabels)):
    ax = axes[i//2, i%2]
    
    # Violin plot
    sns.violinplot(data=df, y=metric, ax=ax, color="#E8E8E8", linewidth=1.2, cut=0)
    
    # Boxplot overlay
    sns.boxplot(data=df, y=metric, ax=ax, width=0.25, 
                boxprops=dict(facecolor="white", edgecolor="#2C3E50", linewidth=1.5),
                medianprops=dict(color="#E74C3C", linewidth=2.5),
                whiskerprops=dict(color="#2C3E50", linewidth=1.5),
                capprops=dict(color="#2C3E50"),
                flierprops=dict(marker='o', color='#2C3E50', alpha=0.6))
    
    # Individual points (stripplot)
    sns.stripplot(data=df, y=metric, ax=ax, color="#2C3E50", 
                  alpha=0.75, jitter=0.18, size=4.5)
plt.savefig("alpha_diversity_plot.png", bbox_inches='tight')
plt.savefig("alpha_diversity_plot.pdf", bbox_inches='tight')

print("✅ Alpha Diversity Figure saved!")