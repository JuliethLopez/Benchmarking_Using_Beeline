import pandas as pd
from pathlib import Path

# List of input GENIE3 output files
esc_dir = Path(
    "/lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs/experimental/ESC"
)

# Search for all GENIE3/outFile.txt files in every subdirectory inside ESC
for method in ["GENIE3","GRNBOOST2"]:
    input_files = list(esc_dir.glob(f"*/{method}/outFile.txt"))

    for file_path in input_files:
        if not file_path.exists():
            print(f"File {file_path} does not exist, skipping...")
            continue

        # Read tab-delimited output file
        df = pd.read_csv(file_path, sep="\t", header=0)

        # Define output CSV path in the same directory as outFile.txt
        output_path = file_path.parent / "rankedEdges.csv"

        # Save with renamed columns
        df_out = df[["TF", "target", "importance"]].rename(
            columns={"TF": "Gene1", "target": "Gene2", "importance": "EdgeWeight"}
        )
        df_out.to_csv(output_path, sep="\t", index=False)

        print(f"Successfully generated: {output_path}")


input_files = list(esc_dir.glob("*/hdWGCNA/outFile.txt"))

for file_path in input_files:
    if not file_path.exists():
        print(f"{file_path} does not exist, skipping...")
        continue
    else:
        print(f"Output file found, proceeding with formatting: {file_path}")

    # Read output
    OutDF = pd.read_csv(file_path, sep="\t", header=0)

    # Format the columns according to hdWGCNA specifications
    required_cols = ["Gene1", "Gene2", "EdgeWeight"]

    if not all(col in OutDF.columns for col in required_cols):
        outFile = OutDF[["tf", "gene", "Cor"]].copy()
        outFile.columns = required_cols
    else:
        outFile = OutDF[required_cols].copy()

    # Save output to rankedEdges.csv in the same directory
    output_path = file_path.parent / "rankedEdges.csv"
    outFile.to_csv(output_path, sep="\t", index=False)

    print(f"Successfully generated: {output_path}")

input_files = list(esc_dir.glob("*/PIDC/outFile.txt"))

for file_path in input_files:
    if not file_path.exists():
        print(f"{file_path} does not exist, skipping...")
        continue

    print(f"Processing PIDC file: {file_path}")

    # Read headerless tab-delimited file
    OutDF = pd.read_csv(file_path, sep=r"\s+", header=None)

    # Select first 3 columns and assign standard column headers
    outFile = OutDF[[0, 1, 2]].copy()
    outFile.columns = ["Gene1", "Gene2", "EdgeWeight"]

    output_path = file_path.parent / "rankedEdges.csv"
    outFile.to_csv(output_path, sep="\t", index=False)

    print(f"Successfully generated: {output_path}")

input_files = list(
    esc_dir.glob(
        "*/INFERELATOR/not_imputation_not_prior/random_seed_0/network.tsv.gz"
    )
)

for file_path in input_files:
    if not file_path.exists():
        print(f"{file_path} does not exist, skipping...")
        continue

    print(f"Processing INFERELATOR file: {file_path}")

    # Read network file
    OutDF = pd.read_csv(file_path, sep="\t", index_col=None)

    # Reformat columns
    outFile = OutDF[["regulator", "target", "model_coefficient"]].copy()
    outFile.columns = ["Gene1", "Gene2", "EdgeWeight"]

    # Write rankedEdges.csv into the base INFERELATOR directory
    inferelator_dir = file_path.parents[2]  # Navigates up to .../INFERELATOR
    output_path = inferelator_dir / "rankedEdges.csv"

    outFile.to_csv(output_path, sep="\t", index=False)
    print(f"Successfully generated: {output_path}")

P_VAL_THRESHOLD = 0.01

# Search for all PPCOR/outFile.txt files
input_files = list(esc_dir.glob("*/PPCOR/outFile.txt"))

for file_path in input_files:
    if not file_path.exists():
        print(f"{file_path} does not exist, skipping...")
        continue

    print(f"Processing PPCOR file: {file_path}")

    # Read tab-delimited output file with header
    OutDF = pd.read_csv(file_path, sep="\t", header=0)

    # Edges with significant p-value
    part1 = OutDF.loc[OutDF["pValue"] <= P_VAL_THRESHOLD].copy()
    part1["absCorVal"] = part1["corVal"].abs()
    part1 = part1.sort_values("absCorVal", ascending=False)

    # Standardize part1 columns
    part1_out = part1[["Gene1", "Gene2", "corVal"]].rename(
        columns={"corVal": "EdgeWeight"}
    )

    # Edges without significant p-value (EdgeWeight set to 0)
    part2 = OutDF.loc[OutDF["pValue"] > P_VAL_THRESHOLD].copy()
    part2_out = part2[["Gene1", "Gene2"]].copy()
    part2_out["EdgeWeight"] = 0

    # Combine both parts
    outFile = pd.concat([part1_out, part2_out], ignore_index=True)

    output_path = file_path.parent / "rankedEdges.csv"
    outFile.to_csv(output_path, sep="\t", index=False)

    print(f"Successfully generated: {output_path}")