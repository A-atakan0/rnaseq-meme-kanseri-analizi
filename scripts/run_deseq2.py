import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Veriyi yukle
counts = pd.read_csv("results/count_matrix.csv", index_col=0).T  # ornek x gen olacak sekilde transpoze
metadata = pd.read_csv("raw_data/metadata.csv", index_col=0)

# Sirlamayi eslestir
metadata = metadata.loc[counts.index]

# Dusuk expresyonlu genleri filtrele (en az 10 toplam okuma)
counts = counts.loc[:, counts.sum(axis=0) >= 10]

print(f"Analiz edilecek gen sayisi: {counts.shape[1]}")
print(f"Ornek sayisi: {counts.shape[0]}")
print(metadata)

# DESeq2 modelini kur
dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
)
dds.deseq2()

# Istatistiksel test: tumor vs normal
stat_res = DeseqStats(dds, contrast=["condition", "tumor", "normal"])
stat_res.summary()

results_df = stat_res.results_df
results_df = results_df.sort_values("padj")
results_df.to_csv("results/deseq2_results.csv")

print("\nEn anlamli 10 gen:")
print(results_df.head(10))

sig_genes = results_df[(results_df["padj"] < 0.05) & (abs(results_df["log2FoldChange"]) > 1)]
print(f"\nAnlamli farkli eksprese olan gen sayisi (padj<0.05, |log2FC|>1): {len(sig_genes)}")
