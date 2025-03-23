import pandas as pd
import glob
import os
import numpy as np

# Base directory containing all subdirectories
base_directory = os.path.join("..", "simulations", "Tsitsi_0,9_3,0_2,0")

# Columns to ignore
ignore_columns = {"Xlabel", "Ylabel", "Radius", "TimeStep"}  # Replace with actual column names

# Iterate over each subdirectory in the base directory
for subdir in sorted(os.listdir(base_directory)):
    subdir_path = os.path.join(base_directory, subdir)

    # Check if it's a directory
    if not os.path.isdir(subdir_path):
        continue

    # File pattern for all the data files in this subdirectory
    file_pattern = os.path.join(subdir_path, "connections_dat_*")
    files = sorted(glob.glob(file_pattern))

    # Check if files are found
    if not files:
        print(f"No files found in directory {subdir_path}, skipping.")
        continue

    # Load column names from the first file
    sample_df = pd.read_csv(files[0], sep='\s+')
    columns = [col for col in sample_df.columns if col not in ignore_columns]  # Filter out ignored columns

    # Create a dictionary to store the timeSteps and data for statistical analysis
    time_steps = sample_df.index.values  # Extract time steps from the index (row numbers)
    column_data = {col: [] for col in columns}

    # Read each file and store data for each column
    for file in files:
        df = pd.read_csv(file, sep='\s+')
        for col in columns:  # Process only required columns
            column_data[col].append(df[col].values)

    # Create an output directory for the current subdirectory
    output_directory = os.path.join(base_directory, f"{subdir}_statistics")
    os.makedirs(output_directory, exist_ok=True)

    # Perform statistical analysis and write results for each column
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
        print(f"Saved statistics for column '{col}' in {subdir} to {output_file}")
