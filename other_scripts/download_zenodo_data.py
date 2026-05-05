# remember to cite zenodo_get
from zenodo_get import download

# Download all files from a record

# (McCalla et al., 2023b) https://zenodo.org/records/5909090
# We will use a subset of Human pluripotent stem cells (hPSCs) or dataset HAN. Those are Emrbionic stem cell (ESC)
download(record_or_doi = "10.5281/zenodo.5909090", output_dir="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/5909090")
#file_glob="*.csv"

# (Li et al., 2025) https://zenodo.org/records/15711227
# We will use subsets of PBMCs. Those are monocytes, DC, B, CD4 T
download(record_or_doi = "10.5281/zenodo.15711227", output_dir="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/15711227")
# The reference is NOT AVAILABLE, it needs to be constructed from CISTROME (ChiP-seq data)
download(record_or_doi = "10.5281/zenodo.15720607", output_dir="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/15720607")

# (Pratapa et al., 2020) https://zenodo.org/records/19009603
# We will use hHep, hESC
download(record_or_doi = "10.5281/zenodo.19009603", output_dir="/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/data/experimental/19009603")
