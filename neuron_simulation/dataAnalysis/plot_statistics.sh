#!/bin/bash

# Directory where statistics files are stored
output_dir="simulations/lattice_150_statistics"

# Loop over each _statistics.txt file in the output directory
for file in $output_dir/*_statistics.txt; do
    # Extract the column name (remove the '_statistics.txt' part)
    column_name=$(basename "$file" _statistics.txt)
    
    # Set the output PNG file name
    output_file="${output_dir}/${column_name}_plot.png"

    # Generate the Gnuplot script dynamically
    gnuplot -e "
        set terminal pngcairo size 1280,720 enhanced font 'Arial,12';
        set output '${output_file}';
        set title 'Statistics for ${column_name}' font ',14';
        set xlabel 'Time Steps' font ',12';
        set ylabel 'Mean ± Std Dev' font ',12';
        set grid;
        set key outside;

        plot '${file}' using 1:2:3 every 5000 with yerrorbars lc rgb 'red' pt 7 ps 1.5 lt 1 lw 1 title 'Mean ± Std Dev', \
             '' using 1:2 with lines lc rgb 'red' lt 1 lw 3 title 'Mean';
    "
    echo "Plot saved as ${output_file}"
done