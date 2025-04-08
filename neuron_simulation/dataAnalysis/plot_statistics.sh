#!/bin/bash

# Base directory where the statistics folders are stored
base_dir="../simulations/Ising_0,9_10,0_6,0"

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
    # Terminal settings
    set terminal pngcairo size 1280,720 enhanced font 'Arial,14'
    set output '${output_file}'
    
    # Color definitions
    set style line 1 lc rgb '#648FFF' lt 1 lw 3 pt 7 ps 0   # IBM Blue (Primary data)
    set style line 2 lc rgb '#785EF0' lt 1 lw 3 pt 7 ps 0   # IBM Purple (Secondary)
    set style line 3 lc rgb '#DC267F' lt 1 lw 3 pt 7 ps 0   # IBM Red (Tertiary)
    set style line 4 lc rgb '#FE6100' lt 1 lw 3 pt 7 ps 0   # IBM Orange (Highlight)
    set style line 5 lc rgb '#009E73' lt 1 lw 3 pt 7 ps 0   # IBM Cyan (Complementary)


    # Layout - adjust margins to give more space for labels
    set lmargin at screen 0.12
    set rmargin at screen 0.95
    set bmargin at screen 0.12
    set tmargin at screen 0.95
    
    # Axis and ticks - complete box with ticks pointing in
    set border linewidth 1.5
    set xtics scale 2, 1
    set ytics scale 2,1
    set tics in  # ticks point inward
    
    # Set logarithmic scale for both axes    
    set logscale x
    set logscale y
    
    # Only use scientific notation when needed
    set autoscale
    stats '${file}' using 2 nooutput
    y_range = STATS_max - STATS_min
    if (y_range < 1000 && y_range > 0) {
        set format y '%g'  # Use normal format for reasonable ranges
    }
    stats '${file}' using 1 nooutput
    x_range = STATS_max - STATS_min
    if (x_range < 1000 && x_range > 0) {
        set format x '%g'  # Use normal format for reasonable ranges
    }
    
    # Title and labels - with more offset from axes
    #set title '${column_name}' font ',16' offset 0,1
    set xlabel '{/Symbol t}' font 'Times New Roman,28' offset 0,-0.2
    set ylabel '${column_name}' font 'Times New Roman,28' offset -0.5,0
    
    # Legend placement
    set key at graph 0.98, graph 0.75 
    set key font 'Times New Roman, 16'
    set key spacing 1.8                 # Increased line spacing (default 1.25)
    set key width +0                       # Extra width allowance (character units)
    set key height +0                   # Extra height per entry (character units)
    set key box lt -1 lw 1              # Stronger border (optional)
    set key maxrows 5                   # Limit rows before column break (optional)
    
    # Confidence interval style
    set style fill transparent solid 0.25 noborder
    
    # Get maximum value from last data point
    stats '${file}' using 2 nooutput
    max_value = STATS_max
    
    # Plot commands
    plot '${file}' using 1:(\$2-\$3):(\$2+\$3) with filledcurves ls 1 notitle, \
         '${file}' using 1:2 with lines ls 1 lw 3 title '${column_name}', \
         max_value with lines ls 2 lw 2 dt 2 title sprintf("Max: %.2g", max_value)
EOF

        # Run GNUPlot using the .gp file
        gnuplot "$gp_file"

        echo "Plot saved as ${output_file}"
        echo "Gnuplot script saved as ${gp_file}"
    done
done