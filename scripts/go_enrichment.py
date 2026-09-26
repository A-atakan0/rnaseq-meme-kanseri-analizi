import pandas as pd
import gseapy as gp

results = pd.read_csv("results/deseq2_results_annotated.csv", index_col=0)
sig_up = results[(results['padj']<0.05) & (results['log2FoldChange']>1)]['gene_symbol'].dropna().tolist()
sig_down = results[(results['padj']<0.05) & (results['log2FoldChange']<-1)]['gene_symbol'].dropna().tolist()

print(f"Yukari regule gen sayisi: {len(sig_up)}")
print(f"Asagi regule gen sayisi: {len(sig_down)}")

enr_up = gp.enrichr(gene_list=sig_up, gene_sets='GO_Biological_Process_2021', outdir=None)
enr_up.results.sort_values('Adjusted P-value').head(10).to_csv("results/GO_upregulated_top10.csv")
print("\nYukari regule genlerde en anlamli 10 GO terimi:")
print(enr_up.results.sort_values('Adjusted P-value')[['Term','Adjusted P-value']].head(10).to_string())

enr_down = gp.enrichr(gene_list=sig_down, gene_sets='GO_Biological_Process_2021', outdir=None)
enr_down.results.sort_values('Adjusted P-value').head(10).to_csv("results/GO_downregulated_top10.csv")
print("\nAsagi regule genlerde en anlamli 10 GO terimi:")
print(enr_down.results.sort_values('Adjusted P-value')[['Term','Adjusted P-value']].head(10).to_string())
