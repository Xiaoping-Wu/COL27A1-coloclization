import pandas as pd

# Snakemake automatically passes input/output variables
input_files = snakemake.input.files
output_file = snakemake.output[0]

# Read and merge (row-wise)
dfs = [pd.read_csv(f, sep="\t") for f in input_files]
merged = pd.concat(dfs, ignore_index=True)

# Apply filtering criteria (example: p < 5e-8)
filtered = merged[(merged["PP.H3"] >= 0.7) | (merged["PP.H4"] >= 0.7)]

# Save result
filtered.to_csv(output_file, sep="\t", index=False)

print(f"Merged {len(input_files)} files → {len(filtered)} rows kept")