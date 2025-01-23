#!/bin/bash

# Directories to process
directories=("simulations/random_100_statistics" "simulations/random_125_statistics" "simulations/random_150_statistics")

# Output directory for the combined plots
output_dir="simulations/random_combined_plots"
mkdir -p "$output_dir"

# Get the list of statistics files from the first directory
for file in "${directories[0]}"/*_statistics.txt; do
    # Extract the base file name (without directory and extension)
    base_name=$(basename "$file" _statistics.txt)
    
    # Set the output PNG file name
    output_file="${output_dir}/${base_name}_combined_plot.png"
    
    # Generate the Gnuplot command dynamically
    gnuplot -e "
        set terminal pngcairo size 1280,720 enhanced font 'Times-Roman,20';
        set output '${output_file}';
        set title '${base_name} vs Time' font 'Times-Roman,20';
        set xlabel 'Time Steps' font 'Times-Roman,20';
        set ylabel '${base_name}' font 'Times-Roman,20';
        set grid;
        set key at graph 1,0.25;
        set key box;
        set key opaque;

        set format x '10^{%L}';

        plot '${directories[0]}/${base_name}_statistics.txt' using 1:2:3 every 5000 with yerrorbars lc rgb 'red' pt 7 ps 1.5 lt 1 lw 1 notitle, \
             '${directories[0]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'red' lt 1 lw 3 title '0.39 {/Symbol F}', \
             '${directories[1]}/${base_name}_statistics.txt' using 1:2:3 every 5000 with yerrorbars lc rgb 'blue' pt 7 ps 1.5 lt 1 lw 1 notitle, \
             '${directories[1]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'blue' lt 1 lw 3 title '0.25 {/Symbol F}', \
             '${directories[2]}/${base_name}_statistics.txt' using 1:2:3 every 5000 with yerrorbars lc rgb 'green' pt 7 ps 1.5 lt 1 lw 1 notitle, \
             '${directories[2]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'green' lt 1 lw 3 title '0.17 {/Symbol F}';
    "
    echo "Combined plot saved as ${output_file}"
done
