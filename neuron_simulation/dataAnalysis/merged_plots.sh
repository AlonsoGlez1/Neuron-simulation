#!/bin/bash

# Base directory where the statistics folders are stored
base_dir="../simulations"

# Define groups to process
groups=("random" "lattice" "hexagonal" "cluster")

# Loop over each group
for group in "${groups[@]}"; do
    echo "Processing group: $group"

    # Find all matching directories for the current group (e.g., random_100_statistics, random_125_statistics, ...)
    directories=($(ls -d "$base_dir/${group}_"*_statistics 2>/dev/null))

    # Skip if no directories are found
    if [ ${#directories[@]} -eq 0 ]; then
        echo "No directories found for $group, skipping..."
        continue
    fi

    # Output directory for the combined plots
    output_dir="${base_dir}/${group}_combined_plots"
    mkdir -p "$output_dir"

    # Get the list of statistics files from the first directory
    for file in "${directories[0]}"/*_statistics.txt; do
        # Extract the base file name (without directory and extension)
        base_name=$(basename "$file" _statistics.txt)

        # Set the output PNG file name
        output_file="${output_dir}/${base_name}_combined_plot.png"

        # Start Gnuplot command
        gnuplot_cmd="
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
            set style fill transparent solid 0.25;
            set style fill noborder;
        "

        # Define colors for each dataset
        colors=("red" "blue" "green" "purple" "orange")

        # Loop over directories and add plotting commands
        plot_cmd=""
        for i in "${!directories[@]}"; do
            dir="${directories[$i]}"
            color="${colors[$i % ${#colors[@]}]}"  # Cycle through colors if more than 5 directories
            phi_index=$((i + 1))  # Phi index starts from 1

            # Append plot lines for this directory
            plot_cmd+="'${dir}/${base_name}_statistics.txt' using 1:(\$2-\$3):(\$2+\$3) with filledcurves lc rgb '${color}' notitle, "
            plot_cmd+="'${dir}/${base_name}_statistics.txt' using 1:2 with lines lc rgb '${color}' lt 1 lw 3 title '{/Symbol F}_${phi_index}', "
        done

        # Remove the last comma
        plot_cmd=${plot_cmd%, }

        # Complete the Gnuplot command
        gnuplot_cmd+="plot ${plot_cmd};"

        # Run Gnuplot
        gnuplot -e "$gnuplot_cmd"

        echo "Combined plot saved as ${output_file}"
    done
done
