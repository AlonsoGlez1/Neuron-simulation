#!/bin/bash

# Base directory where the statistics folders are stored
base_dir="../simulations/Ising_0,9_3,0_3,0"

# Find all directories matching the pattern *_statistics
for output_dir in "$base_dir"/*_statistics; do
    # Ensure it's a directory
    [ -d "$output_dir" ] || continue 

    echo "Processing directory: $output_dir"

    # Loop over each _statistics.txt file in the directory
    for file in "$output_dir"/*_statistics.txt; do
        # Ensure the file exists (avoid processing wildcard if no match)
        [ -f "$file" ] || { echo "Skipping non-existent file: $file"; continue; }

        # Extract the column name
        column_name=$(basename "$file" _statistics.txt)

        # Set the output PNG file name
        output_file="${output_dir}/${column_name}_plot.png"

        # Set the output .gp file name
        gp_file="${output_dir}/${column_name}_plot.gp"

        # Debugging: print values
        echo "Processing file: $file"
        echo "Output plot: $output_file"
        echo "Output .gp file: $gp_file"

        # Create the .gp file with the Gnuplot commands
        cat <<EOF > "$gp_file"
set terminal pngcairo size 1280,720 enhanced font 'Arial,12';
set output '${output_file}';
set title 'Statistics for ${column_name}' font ',14';
set xlabel 'Time Steps' font ',12';
set ylabel 'Mean ± Std Dev' font ',12';
set grid;
set key outside;
set style fill transparent solid 0.25; 
set style fill noborder;

plot '${file}' using 1:(\$2-\$3):(\$2+\$3) with filledcurves lc 'blue' notitle, \
     '${file}' using 1:2 with lines lc 'blue' lt 1 lw 3 title 'Mean';
EOF

        # Run GNUPlot using the .gp file
        gnuplot "$gp_file"

        echo "Plot saved as ${output_file}"
        echo "Gnuplot script saved as ${gp_file}"
    done
done