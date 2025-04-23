import pandas as pd

# Step 1: Read the file into a DataFrame (column 11 is loaded as string to avoid type issues)
file_path = 'scan.dat.gz'
nf = pd.read_csv(file_path, sep=r'\s+', compression='gzip', dtype={11: str})

# Step 2: Identify the last value of column 11
last_value = nf.iloc[-1, 11]

# Step 3: Remove the 'ls' from the last value (if present)
if isinstance(last_value, str) and last_value.endswith('ls'):
    nf.iloc[-1, 11] = last_value[:-2]  # Remove the last two characters ('ls')

# Step 4: Save the updated DataFrame back to a file (new file to preserve the original)
output_file_path = 'scan_fixed.dat.gz'
nf.to_csv(output_file_path, sep='\t', compression='gzip', index=False)

print(f"File saved successfully to {output_file_path}")
