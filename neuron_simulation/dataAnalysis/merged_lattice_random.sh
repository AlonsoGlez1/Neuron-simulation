#!/bin/bash

# Directories to process
directories=("simulations/random_100_statistics" "simulations/lattice_100_statistics" "simulations/random_125_statistics" "simulations/lattice_125_statistics" "simulations/random_150_statistics" "simulations/lattice_150_statistics")

# Output directory for the combined plots
output_dir="simulations/random_combined_random_lattice"
mkdir -p "$output_dir"

# Get the list of statistics files from the first directory
for file in "${directories[0]}"/*_statistics.txt; do
    # Extract the base file name (without directory and extension)
    base_name=$(basename "$file" _statistics.txt)
    
    # Set the output PNG file name
    output_file="${output_dir}/${base_name}_combined_plot.png"
    
    # Generate the Gnuplot command dynamically
    gnuplot -e "
        set terminal pngcairo size 1280,720 enhanced font 'Times-Roman,15';
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
             '${directories[0]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'red' lt 1 lw 3 title '0.39 {/Symbol F} random', \
             '${directories[1]}/${base_name}_statistics.txt' using 1:2:3 every 5000 with yerrorbars lc rgb 'blue' pt 7 ps 1.5 lt 1 lw 1 notitle, \
             '${directories[1]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'blue' lt 1 lw 3 title '0.39 {/Symbol F} lattice', \
             '${directories[4]}/${base_name}_statistics.txt' using 1:2:3 every 5000 with yerrorbars lc rgb 'violet' pt 7 ps 1.5 lt 1 lw 1 notitle, \
             '${directories[4]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'violet' lt 1 lw 3 title '0.17 {/Symbol F} random', \
             '${directories[5]}/${base_name}_statistics.txt' using 1:2:3 every 5000 with yerrorbars lc rgb 'aquamarine' pt 7 ps 1.5 lt 1 lw 1 notitle, \
             '${directories[5]}/${base_name}_statistics.txt' using 1:2 with lines lc rgb 'aquamarine' lt 1 lw 3 title '0.17 {/Symbol F} lattice';
    "
    echo "Combined plot saved as ${output_file}"
done
