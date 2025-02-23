import pandas as pd

# Load the CSV files
meta = pd.read_csv(r"C:\Users\fonte\OneDrive - UGent\AJ 2024-2025\Master's_Dissertation\AML-CITEseq\Data\metadata_after_merge.csv",
                    header=0, sep=",", index_col=0)

patientoverview = pd.read_csv(r"C:\Users\fonte\OneDrive - UGent\AJ 2024-2025\Master's_Dissertation\AML-CITEseq\Data\Metadata_table_clean.csv",
                              header=0, sep=";")

# Extract SABxxx from the 'orig.ident' column in meta
meta['exp'] = meta['orig.ident'].str.extract(r'(SAB\d{3})')

# Merge 'meta' with 'patientoverview' on 'exp' (assuming 'exp' corresponds to a column in patientoverview)
merged_meta = meta.merge(patientoverview[['exp', 'day', 'Sorting']], on='exp', how='left')

# Print to verify
print(merged_meta[['orig.ident', 'exp', 'day', 'Sorting']].head())  # Display first few rows

# Save the updated meta dataframe (if needed)
merged_meta.to_csv(r"C:\Users\fonte\OneDrive - UGent\AJ 2024-2025\Master's_Dissertation\AML-CITEseq\Data\updated_metadata.csv", index=False)
