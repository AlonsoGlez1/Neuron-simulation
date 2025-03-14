import pandas as pd
import glob
import os
import numpy as np

# Path to directory containing the files
file_directory = os.path.join("..", "simulations", "5k_cluster_150Lx_0,5rad_0,17453pac_10e5time_1branch_2synapses") # Relative path

# File pattern for all the data files
file_pattern = os.path.join(file_directory, "connections_dat_*")
files = sorted(glob.glob(file_pattern))

# Check if files are found
if not files:
    raise FileNotFoundError(f"No files found in directory {file_directory} matching pattern {file_pattern}")

# Load column names from the first file
sample_df = pd.read_csv(files[0], sep = '\s+')
columns = sample_df.columns

# Create a dictionary to store the timeSteps and data for statistical analysis
time_steps = sample_df.index.values  # Extract time steps from the index (row numbers)
column_data = {col: [] for col in columns}

# Read each file and store data for each column
for file in files:
    df = pd.read_csv(file, sep = '\s+')
    for col in columns:
        column_data[col].append(df[col].values)

# Perform statistical analysis and write results for each column
output_directory = os.path.join("..", "simulations", "cluster_150_statistics")
os.makedirs(output_directory, exist_ok=True)

for col, data in column_data.items():
    # Combine all rows for this column across files
    combined_data = np.column_stack(data)  # Each row corresponds to a time step, columns to files
    
    # Calculate mean and standard deviation for each time step
    mean = np.mean(combined_data, axis=1)
    std_dev = np.std(combined_data, axis=1)
    
    # Calculate Standard Error of the Mean (SEM)
    N = combined_data.shape[1]  # Number of simulations (files)
    sem = std_dev / np.sqrt(N)
    
    # Create a DataFrame for output
    results_df = pd.DataFrame({
        "timeSteps": time_steps,
        "mean": mean,
        "std_dev": std_dev,
        "sem": sem
    })
    
    # Save to a file named after the column
    output_file = os.path.join(output_directory, f"{col}_statistics.txt")
    results_df.to_csv(output_file, index=False, sep="\t")
    print(f"Saved statistics for column '{col}' to {output_file}")
