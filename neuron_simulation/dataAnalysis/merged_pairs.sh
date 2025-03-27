#!/bin/bash

# Base directory where the statistics folders are stored
base_dir="../simulations/Tsitsi_0,9_3,0_2,0"

# Define the group pairs to compare
comparisons=(
    "random lattice"
    "random cluster"
    "lattice hexagonal"
)

# Loop over each comparison pair
for pair in "${comparisons[@]}"; do
    # Extract group names
    group1=$(echo $pair | awk '{print $1}')
    group2=$(echo $pair | awk '{print $2}')
    
    echo "Processing comparison: $group1 vs $group2"

    # Find all matching directories for each group
    directories1=($(ls -d "$base_dir/"*"${group1}_"*_statistics 2>/dev/null))
    directories2=($(ls -d "$base_dir/"*"${group2}_"*_statistics 2>/dev/null))

    # Skip if either group has no directories
    if [ ${#directories1[@]} -eq 0 ] || [ ${#directories2[@]} -eq 0 ]; then
        echo "Skipping $group1 vs $group2 (missing data)"
        continue
    fi

    # Output directory for the combined plots
    output_dir="${base_dir}/${group1}_vs_${group2}_plots"
    mkdir -p "$output_dir"

    # Get the list of statistics files from the first directory of group1
    for file in "${directories1[0]}"/*_statistics.txt; do
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

            set format x '%.1te%T';
            set style fill transparent solid 0.25;
            set style fill noborder;
        "

        # Define colors for each dataset
        colors=("red" "blue" "green" "purple" "orange" "cyan")

        # Loop over directories and add plotting commands
        plot_cmd=""
        index=1
        for i in "${!directories1[@]}"; do
            dir1="${directories1[$i]}"
            dir2="${directories2[$i]}"
            color1="${colors[$((2 * i % ${#colors[@]}))]}"
            color2="${colors[$(((2 * i + 1) % ${#colors[@]}))]}"

            # Append plot lines for this pair
            plot_cmd+="'${dir1}/${base_name}_statistics.txt' using 1:(\$2-\$3):(\$2+\$3) with filledcurves lc rgb '${color1}' notitle, "
            plot_cmd+="'${dir1}/${base_name}_statistics.txt' using 1:2 with lines lc rgb '${color1}' lt 1 lw 3 title '{/Symbol F}_${index} ${group1}', "
            
            plot_cmd+="'${dir2}/${base_name}_statistics.txt' using 1:(\$2-\$3):(\$2+\$3) with filledcurves lc rgb '${color2}' notitle, "
            plot_cmd+="'${dir2}/${base_name}_statistics.txt' using 1:2 with lines lc rgb '${color2}' lt 1 lw 3 title '{/Symbol F}_${index} ${group2}', "

            index=$((index + 1))
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
