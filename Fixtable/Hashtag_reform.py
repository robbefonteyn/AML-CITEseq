import pandas as pd

# Load the CSV file
file_path = "Metadata_table_clean.csv"
df = pd.read_csv(file_path, delimiter=';')

# Replace values in the 'hashtag' column, keeping 'NA' as is
df['hashtag'] = df['hashtag'].astype(str).apply(lambda x: f'Hashtag{x}' if x.strip().upper() != 'NA' else x)

# Save the updated file
output_path = "Metadata_table_clean.csv"  # Update with your desired output file path
df.to_csv(output_path, index=False, sep=';')