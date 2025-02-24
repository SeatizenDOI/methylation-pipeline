import os
import csv
import glob
import sys

# Define the keys we want to extract
keys = [
    "Total number of C's analysed",
    "Total methylated C's in CpG context",
    "Total methylated C's in CHG context",
    "Total methylated C's in CHH context",
    "Total methylated C's in Unknown context",
    "Total unmethylated C's in CpG context",
    "Total unmethylated C's in CHG context",
    "Total unmethylated C's in CHH context",
    "Total unmethylated C's in Unknown context",
    "C methylated in CpG context",
    "C methylated in CHG context",
    "C methylated in CHH context",
    "C methylated in Unknown context (CN or CHN)"
]

# Get directory from command line arguments or use current directory
directory = sys.argv[1] if len(sys.argv) > 1 else "."

# Dictionary to store extracted data
data = {key: {} for key in keys}

# Get all txt files in specified directory
txt_files = glob.glob(os.path.join(directory, "*.txt"))

def trim_filename(filename):
    """Trim filename to remove suffix after 'R1' or 'R2'"""
    parts = filename.split("_R")
    return parts[0] + "_R" + parts[1][0] if len(parts) > 1 else filename

for file in txt_files:
    file_name = os.path.splitext(os.path.basename(file))[0]  # Get file name without extension
    trimmed_name = trim_filename(file_name)
    with open(file, "r") as f:
        for line in f:
            for key in keys:
                if line.startswith(key):
                    value = line.split("\t")[-1].strip()
                    data[key][trimmed_name] = value

# Write to CSV
csv_file = os.path.join(directory, "output.csv")
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    header = ["Metric"] + [trim_filename(os.path.splitext(os.path.basename(f))[0]) for f in txt_files]
    writer.writerow(header)
    
    for key in keys:
        row = [key] + [data[key].get(f, "") for f in header[1:]]
        writer.writerow(row)

print(f"Data successfully written to {csv_file}")
