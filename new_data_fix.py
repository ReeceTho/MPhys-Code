import gzip

def format_to_sci_notation(value):
    try:
        return "{:.6e}".format(float(value))
    except ValueError:
        return value  # Keep non-numeric values (e.g., headers) as-is

def process_and_compress(input_gz_file, output_gz_file):
    with gzip.open(input_gz_file, 'rt') as infile, gzip.open(output_gz_file, 'wt') as outfile:
        for line in infile:
            parts = line.strip().split()
            formatted_parts = [format_to_sci_notation(p) for p in parts]
            outfile.write(" ".join(formatted_parts) + "\n")

# Example usage:
input_gz = "scan.dat.gz"
output_gz = "scan_formatted.dat.gz"
process_and_compress(input_gz, output_gz)
print("done")