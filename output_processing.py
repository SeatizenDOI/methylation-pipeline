import os
import pandas as pd
import argparse

def process_bedgraph_files(input_folder, threshold=0.5):
    """
    Reads bedGraph files from the input folder and constructs a methylation matrix.
    
    Args:
        input_folder (str): Path to the folder containing bedGraph files.
        threshold (float): Minimum percentage of samples a (chr, start) must appear in to be included.
    
    Returns:
        pd.DataFrame: Methylation matrix.
    """
    all_data = {}
    positions = {}

    # Read each bedGraph file and store data
    for file in os.listdir(input_folder):
        if file.endswith(".bedGraph.gz"):
            sample_name = file.split("_QC")[0]  # Use filename as sample name
            file_path = os.path.join(input_folder, file)
            
            df = pd.read_csv(file_path, sep="\t", header=0, names=["Chromosome", "Start", "End", "Methylation"])
            df["Identifier"] = df["Chromosome"] + "_" + df["Start"].astype(str)
            
            print(df)
            # Store sample data
            all_data[sample_name] = df.set_index("Identifier")["Methylation"].to_dict()
            
            # Count occurrences of each identifier
            for identifier in df["Identifier"]:
                positions[identifier] = positions.get(identifier, 0) + 1

    num_samples = len(all_data)
    min_samples = int(threshold * num_samples)  # Calculate threshold count

    # Filter positions based on threshold
    filtered_positions = {pos for pos, count in positions.items() if count >= min_samples}

    # Create a dataframe with selected positions
    methylation_matrix = pd.DataFrame.from_dict(all_data, orient="index", dtype=float)
    methylation_matrix = methylation_matrix.reindex(columns=sorted(filtered_positions)).fillna(0)

    # Reset index and rename the first column to "SampleID"
    methylation_matrix.reset_index(inplace=True)
    methylation_matrix.rename(columns={"index": "SampleID"}, inplace=True)

    return methylation_matrix

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process bedGraph files into a methylation matrix.")
    parser.add_argument("input_folder", type=str, help="Folder containing bedGraph files")
    parser.add_argument("--threshold", type=float, default=0.5, help="Minimum percentage of samples a position must appear in")

    args = parser.parse_args()
    
    matrix = process_bedgraph_files(args.input_folder, args.threshold)
    output_path = os.path.join(args.input_folder, "methylation_matrix.csv")
    matrix.to_csv(output_path, index=False, sep=",")

    print(f"Methylation matrix saved to {output_path}")
