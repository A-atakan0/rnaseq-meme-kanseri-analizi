import pandas as pd
import glob

files = sorted(glob.glob("results/*_counts.txt"))
count_data = {}

for f in files:
    sample = f.split("/")[-1].replace("_counts.txt", "")
    df = pd.read_csv(f, sep="\t", comment="#")
    # featureCounts ciktisinda son sutun okuma sayisi, ilk sutun Geneid
    df = df[["Geneid", df.columns[-1]]]
    df.columns = ["Geneid", sample]
    count_data[sample] = df.set_index("Geneid")[sample]

matrix = pd.DataFrame(count_data)
matrix.to_csv("results/count_matrix.csv")
print("Sayim matrisi olusturuldu:", matrix.shape, "gen x ornek")
print(matrix.head())
