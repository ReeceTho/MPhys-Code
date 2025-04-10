import pandas as pd
import shutil

# File paths
original_file = 'scan.dat.gz'
backup_file = 'scan_backup.dat.gz'

# Backup original
shutil.copy(original_file, backup_file)

# Read from gzip with header detection
df = pd.read_csv(original_file, sep='\s+', compression='gzip', low_memory=False)

# Convert MD1 to float safely
df['MD1'] = pd.to_numeric(df['MD1'], errors='coerce')
df = df.dropna(subset=['MD1'])
df = df[df['MD1'] >= 10]

# Write back to gzip
df.to_csv(original_file, index=False, sep=' ', compression='gzip')
print("finished")
