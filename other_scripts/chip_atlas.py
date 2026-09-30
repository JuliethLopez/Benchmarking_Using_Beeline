import requests
import pandas as pd
from pathlib import Path
import os

os.chdir("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data")

####################################################
# Get TFs for one cell type
####################################################

def get_tfs(cell_type=None, genome="hg38"):
    url = "https://chip-atlas.org/data/chip_antigen"
    params = {
        "genome": genome,
        "agClass": "TFs and others"
    }
    if cell_type is not None:
        params["clClass"] = cell_type
    r = requests.get(url, params=params)
    r.raise_for_status()
    return r.json()

####################################################
# Download target genes
####################################################

def download_tf_targets(tf,
                        genome="hg38",
                        window=5,
                        outdir="targets"):
    Path(outdir).mkdir(exist_ok=True)
    url = f"https://chip-atlas.dbcls.jp/data/{genome}/target/{tf}.{window}.tsv"
    outfile = Path(outdir) / f"{tf}.tsv"
    r = requests.get(url)
    if r.status_code == 200:
        outfile.write_bytes(r.content)
        return outfile
    return None

####################################################
# Merge all TF-target pairs
####################################################

def merge_targets(target_folder):
    dfs = []
    for file in Path(target_folder).glob("*.tsv"):
        tf = file.stem
        try:
            df = pd.read_csv(file, sep="\t")
            gene_col = df.columns[0]
            tmp = pd.DataFrame({
                "TF": tf,
                "Target": df[gene_col]
            })
            dfs.append(tmp)
        except Exception as e:
            print(file, e)
    return pd.concat(dfs, ignore_index=True)

####################################################
# Get TFs
####################################################

# All TFs #https://dbarchive.biosciencedbc.jp/en/chip-atlas/data-4.html
zip_file = "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/julieth/chip_atlas_antigen_list.zip"
df = pd.read_csv(
    zip_file,
    sep=",",
    compression="zip"
)

print(df.head())
print(df.columns)

human_tfs = (
    df[
        (df["Genome assembly"] == "hg38") &
        (df["Antigen class"] == "TFs and others")
    ]["Antigen"]
    .drop_duplicates()
    .sort_values()
    .tolist()
)

print(len(human_tfs))
print(human_tfs[:20])

# Beeline table is bigger, so I will continue working with it
df_network = pd.read_csv("/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603/Networks/human/Non-specific-ChIP-seq-network.csv")
print(df_network.columns)
print(len(df_network["Gene1"].unique()))

# Get TF cell type specific
liver_tfs = get_tfs("Liver")
print(liver_tfs)

esc_tfs = get_tfs("Embryo")
print(esc_tfs)

####################################################
# Download all target files
####################################################

for i in range(len(json.loads(esc_tfs))):
    tf = json.loads(esc_tfs)[i]["id"]
    print(f"Downloading targets for {tf}")
    download_tf_targets(tf)

merged = merge_targets("targets")

merged.to_csv("esc_tf_targets.tsv",
              sep="\t",
              index=False)
