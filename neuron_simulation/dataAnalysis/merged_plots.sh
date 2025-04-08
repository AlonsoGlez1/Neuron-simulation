#!/bin/bash

# Base directory where the statistics folders are stored
base_dir="../simulations/Ising_0,9_10,0_6,0"

# Define groups to process
groups=("random" "lattice" "hexagonal" "cluster")

# Loop over each group
for group in "${groups[@]}"; do
    echo "Processing group: $group"

    # Find all matching directories for the current group
    directories=($(ls -d "$base_dir"/*"${group}_"*_statistics 2>/dev/null))

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
        # Extract the base file name
        base_name=$(basename "$file" _statistics.txt)

        # Set the output file names
        output_file="${output_dir}/${base_name}_combined_plot.png"
        gp_file="${output_dir}/${base_name}_combined_plot.gp"

        # Create the .gp file header
cat <<EOF > "$gp_file"
    # Terminal settings
    set terminal pngcairo size 1280,720 enhanced font 'Arial,14'
    set output '${output_file}'

    # Color definitions (IBM color palette)
    set style line 1 lc rgb '#648FFF' lt 1 lw 3 pt 7 ps 0   # IBM Blue
    set style line 2 lc rgb '#785EF0' lt 1 lw 3 pt 7 ps 0   # IBM Purple
    set style line 3 lc rgb '#DC267F' lt 1 lw 3 pt 7 ps 0   # IBM Red
    set style line 4 lc rgb '#FE6100' lt 1 lw 3 pt 7 ps 0   # IBM Orange
    set style line 5 lc rgb '#009E73' lt 1 lw 3 pt 7 ps 0   # IBM Cyan

    # Layout
    set lmargin at screen 0.12
    set rmargin at screen 0.95
    set bmargin at screen 0.12
    set tmargin at screen 0.95

    # Axis and ticks - complete box with ticks pointing in
    set border linewidth 1.5
    set xtics scale 2, 1
    set ytics scale 2, 1
    set tics in  # ticks point inward

    # Set logarithmic scale for both axes    
    set logscale x
    set logscale y

    # Format handling
    set autoscale
    stats '${file}' using 2 nooutput
    y_range = STATS_max - STATS_min
    if (y_range < 1000 && y_range > 0) { set format y '%g' }
    stats '${file}' using 1 nooutput
    x_range = STATS_max - STATS_min
    if (x_range < 1000 && x_range > 0) { set format x '%g' }

    # Labels
    set xlabel '{/Symbol t}' font 'Times New Roman,28' offset 0,-0.2
    set ylabel '${base_name}' font 'Times New Roman,28' offset -0.5,0

    # Legend placement
    set key at graph 0.98, graph 0.75 
    set key font 'Times New Roman, 16'
    set key spacing 1.25                 # Increased line spacing (default 1.25)
    set key width +3                       # Extra width allowance (character units)
    set key height +0                   # Extra height per entry (character units)
    set key box lt -1 lw 1              # Stronger border (optional)
    set key maxrows 5                   # Limit rows before column break (optional)

    # Confidence interval style
    set style fill transparent solid 0.25 noborder

    # Plot command
    plot \\
EOF

        # Add plot entries for each directory
        for i in "${!directories[@]}"; do
            dir="${directories[$i]}"
            color_index=$((i % 5 + 1))
            phi_index=$((i + 1))
            
            # For all but the first entry, add a comma and continuation
            [ $i -gt 0 ] && echo ",\\" >> "$gp_file"
            
            # Add the plot entry
            cat <<EOF >> "$gp_file"
    '${dir}/${base_name}_statistics.txt' using 1:(\$2-\$3):(\$2+\$3) with filledcurves ls ${color_index} notitle, \\
    '' using 1:2 with lines ls ${color_index} lw 3 title '{/Symbol F}_${phi_index}' \\
EOF
        done

        # Remove the trailing backslash from the last line
        sed -i '$ s/ \\$//' "$gp_file"

        # Run Gnuplot
        gnuplot "$gp_file"

        echo "Combined plot saved as ${output_file}"
        echo "Gnuplot script saved as ${gp_file}"
    done
done