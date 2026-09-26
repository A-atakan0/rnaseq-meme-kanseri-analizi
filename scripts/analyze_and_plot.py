import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import mygene

# ---- Sonuclari yukle ----
results = pd.read_csv("results/deseq2_results.csv", index_col=0)
counts = pd.read_csv("results/count_matrix.csv", index_col=0)
metadata = pd.read_csv("raw_data/metadata.csv", index_col=0)

# ---- Gen ID -> Gen ismi cevirisi ----
mg = mygene.MyGeneInfo()
ensembl_ids = results.index.str.split('.').str[0].tolist()  # versiyon numarasini kaldir
gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
id_to_symbol = {g['query']: g.get('symbol', g['query']) for g in gene_info}
results['gene_symbol'] = [id_to_symbol.get(i.split('.')[0], i) for i in results.index]
results.to_csv("results/deseq2_results_annotated.csv")

print("En anlamli 15 gen (isimleriyle):")
print(results[['gene_symbol','log2FoldChange','padj']].head(15))

# ---- VOLCANO PLOT ----
plt.figure(figsize=(9,7))
results['neglog10padj'] = -np.log10(results['padj'].replace(0, 1e-300))
colors = np.where((results['padj']<0.05) & (results['log2FoldChange']>1), '#e74c3c',
          np.where((results['padj']<0.05) & (results['log2FoldChange']<-1), '#3498db', '#bdc3c7'))
plt.scatter(results['log2FoldChange'], results['neglog10padj'], c=colors, s=8, alpha=0.6)
plt.axvline(1, color='gray', linestyle='--', linewidth=0.8)
plt.axvline(-1, color='gray', linestyle='--', linewidth=0.8)
plt.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8)
plt.xlabel("log2 Fold Change (Tumor / Normal)")
plt.ylabel("-log10(padj)")
plt.title("Volcano Plot: Tumor vs Normal Meme Dokusu")
plt.tight_layout()
plt.savefig("figures/volcano_plot.png", dpi=300)
plt.close()
print("Volcano plot kaydedildi.")

# ---- PCA PLOT ----
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

log_counts = np.log2(counts.T + 1)
scaled = StandardScaler().fit_transform(log_counts)
pca = PCA(n_components=2)
pcs = pca.fit_transform(scaled)

pca_df = pd.DataFrame(pcs, columns=['PC1','PC2'], index=counts.columns)
pca_df['condition'] = metadata.loc[pca_df.index, 'condition']

plt.figure(figsize=(7,6))
for cond, color in [('tumor','#e74c3c'), ('normal','#3498db')]:
    sub = pca_df[pca_df['condition']==cond]
    plt.scatter(sub['PC1'], sub['PC2'], label=cond, s=100, color=color, edgecolor='black')
for idx, row in pca_df.iterrows():
    plt.annotate(idx.replace('SRR667',''), (row['PC1'], row['PC2']), fontsize=7)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA: Ornekler Arasi Genel Varyasyon")
plt.legend()
plt.tight_layout()
plt.savefig("figures/pca_plot.png", dpi=300)
plt.close()
print("PCA plot kaydedildi.")

# ---- HEATMAP (en anlamli 30 gen) ----
top30 = results.dropna(subset=['padj']).sort_values('padj').head(30)
heat_data = np.log2(counts.loc[top30.index] + 1)
heat_data.index = top30['gene_symbol'].values
heat_data = heat_data[metadata.sort_values('condition').index]

plt.figure(figsize=(10,10))
col_colors = metadata.loc[heat_data.columns, 'condition'].map({'tumor':'#e74c3c','normal':'#3498db'})
sns.clustermap(heat_data, cmap='vlag', z_score=0, col_colors=col_colors,
               col_cluster=False, figsize=(10,10), yticklabels=True)
plt.savefig("figures/heatmap_top30.png", dpi=300)
plt.close()
print("Heatmap kaydedildi.")
