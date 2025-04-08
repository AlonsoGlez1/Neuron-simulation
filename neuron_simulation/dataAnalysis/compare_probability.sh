#!/bin/bash

# Define the main directories
ISING_DIR="../simulations/Ising_0,9_3,0_6,0"
TSITSI_DIR1="../simulations/Tsitsi_0,9_3,0_1,0"
TSITSI_DIR2="../simulations/Tsitsi_0,9_3,0_0,1"  
OUTPUT_DIR="../simulations/comparison_plots_gamma_6,0_1,0_0,1"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over subdirectories in ISING_DIR that end with "_statistics"
for ISING_SUBDIR in "$ISING_DIR"/*_statistics/; do
    SUBDIR_NAME=$(basename "$ISING_SUBDIR")
    
    # Extract the dynamic part of the subdirectory name
    SHORT_NAME=$(echo "$SUBDIR_NAME" | sed -E 's/.*(cluster|lattice|hexagonal|random)_([^_]+).*/\1_\2/')
    
    TSITSI_SUBDIR1="$TSITSI_DIR1/$SUBDIR_NAME"
    TSITSI_SUBDIR2="$TSITSI_DIR2/$SUBDIR_NAME"
    
    # Create a separate directory for the origin
    ORIGIN_DIR="$OUTPUT_DIR/$SHORT_NAME"
    mkdir -p "$ORIGIN_DIR"
    
    # Check if both matching subdirectories exist
    if [ -d "$TSITSI_SUBDIR1" ] && [ -d "$TSITSI_SUBDIR2" ]; then
        echo "Comparing: $SUBDIR_NAME"
        
        # Loop over text files in the ISING subdirectory
        for ISING_FILE in "$ISING_SUBDIR"/*_statistics.txt; do
            FILE_NAME=$(basename "$ISING_FILE")
            TSITSI_FILE1="$TSITSI_SUBDIR1/$FILE_NAME"
            TSITSI_FILE2="$TSITSI_SUBDIR2/$FILE_NAME"
            
            # Check if matching files exist in both Tsitsi directories
            if [ -f "$TSITSI_FILE1" ] && [ -f "$TSITSI_FILE2" ]; then
                BASE_NAME=$(basename "$FILE_NAME" _statistics.txt)
                PLOT_FILE="$ORIGIN_DIR/${BASE_NAME}_plot.png"
                GNUPLOT_SCRIPT="$ORIGIN_DIR/${BASE_NAME}_plot.gp"
                
                # Create Gnuplot script with consistent formatting
                cat <<EOF > "$GNUPLOT_SCRIPT"
    # Terminal settings
    set terminal pngcairo size 1280,720 enhanced font 'Arial,14'
    set output '$PLOT_FILE'

    # Color definitions (IBM color palette)
    set style line 1 lc rgb '#648FFF' lt 1 lw 3 pt 7 ps 0   # IBM Blue (6.0)
    set style line 2 lc rgb '#DC267F' lt 1 lw 3 pt 7 ps 0   # IBM Red (1.0)
    set style line 3 lc rgb '#FE6100' lt 1 lw 3 pt 7 ps 0   # IBM Orange (0.1)

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
    stats '$ISING_FILE' using 2 nooutput
    y_range = STATS_max - STATS_min
    if (y_range < 1000 && y_range > 0) { set format y '%g' }
    stats '$ISING_FILE' using 1 nooutput
    x_range = STATS_max - STATS_min
    if (x_range < 1000 && x_range > 0) { set format x '%g' }

    # Labels
    set xlabel '{/Symbol t}' font 'Times New Roman,28' offset 0,-0.2
    set ylabel '$BASE_NAME' font 'Times New Roman,28' offset -0.5,0

    # Legend placement
    set key at graph 0.98, graph 0.75 
    set key font 'Times New Roman, 16'
    set key spacing 1.25
    set key width +3
    set key height +0
    set key box lt -1 lw 1
    set key maxrows 5

    # Confidence interval style
    set style fill transparent solid 0.25 noborder

    # Plot command
    plot '$ISING_FILE' using 1:(\$2-\$3):(\$2+\$3) with filledcurves ls 1 notitle, \\
         '$ISING_FILE' using 1:2 with lines ls 1 lw 3 title 'γ = 6.0', \\
         '$TSITSI_FILE1' using 1:(\$2-\$3):(\$2+\$3) with filledcurves ls 2 notitle, \\
         '$TSITSI_FILE1' using 1:2 with lines ls 2 lw 3 title 'γ = 1.0', \\
         '$TSITSI_FILE2' using 1:(\$2-\$3):(\$2+\$3) with filledcurves ls 3 notitle, \\
         '$TSITSI_FILE2' using 1:2 with lines ls 3 lw 3 title 'γ = 0.1'
EOF
                
                # Run Gnuplot
                gnuplot "$GNUPLOT_SCRIPT"
                echo "Generated plot: $PLOT_FILE"
                echo "Gnuplot script saved as: $GNUPLOT_SCRIPT"
            fi
        done
    fi
done

echo "Comparison complete. Check the '$OUTPUT_DIR' directory for plots."