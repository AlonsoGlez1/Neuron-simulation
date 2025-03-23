#!/bin/bash

# Define the main directories
ISING_DIR="../simulations/Ising_0,9_3,0_3,0"
TSITSI_DIR="../simulations/Tsitsi_0,9_3,0_2,0"
OUTPUT_DIR="../simulations/comparison_plots"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over subdirectories in ISING_DIR that end with "_statistics"
for ISING_SUBDIR in "$ISING_DIR"/*_statistics/; do
    SUBDIR_NAME=$(basename "$ISING_SUBDIR")
    
    # Extract the dynamic part of the subdirectory name (e.g., "cluster_100Lx", "lattice_50Lx", etc.)
    SHORT_NAME=$(echo "$SUBDIR_NAME" | sed -E 's/.*(cluster|lattice|hexagonal|random)_([^_]+).*/\1_\2/')
    
    TSITSI_SUBDIR="$TSITSI_DIR/$SUBDIR_NAME"
    
    # Create a separate directory for the origin (e.g., lattice_150Lx)
    ORIGIN_DIR="$OUTPUT_DIR/$SHORT_NAME"
    mkdir -p "$ORIGIN_DIR"
    
    # Check if the matching subdirectory exists in TSITSI_DIR
    if [ -d "$TSITSI_SUBDIR" ]; then
        echo "Comparing: $SUBDIR_NAME"
        
        # Loop over text files in the ISING subdirectory
        for ISING_FILE in "$ISING_SUBDIR"/*.txt; do
            FILE_NAME=$(basename "$ISING_FILE")
            TSITSI_FILE="$TSITSI_SUBDIR/$FILE_NAME"
            
            # Check if the matching file exists in TSITSI_DIR
            if [ -f "$TSITSI_FILE" ]; then
                PLOT_FILE="$ORIGIN_DIR/${FILE_NAME%.txt}.png"
                GNUPLOT_SCRIPT="$ORIGIN_DIR/plot_${FILE_NAME%.txt}.gp"
                
                # Create Gnuplot script with error bars
                PLOT_TITLE="${FILE_NAME%.txt} vs Time"  # Modify the title dynamically based on the file name
                cat <<EOF > "$GNUPLOT_SCRIPT"
set terminal pngcairo enhanced size 1280,720
set output "$PLOT_FILE"
set title "$PLOT_TITLE" font ",20"
set xlabel "Time Steps" font ",20"
set ylabel "Mean ± Std Dev" font ",20"
set grid
set key outside
set style fill transparent solid 0.25
set style fill noborder

# Plot Ising data with error bars
plot "$ISING_FILE" using 1:(\$2-\$3):(\$2+\$3) with filledcurves lc rgb "blue" notitle, \
     "$ISING_FILE" using 1:2 with lines lc rgb "blue" lt 1 lw 3 title "Ising Mean", \
     "$TSITSI_FILE" using 1:(\$2-\$3):(\$2+\$3) with filledcurves lc rgb "red" notitle, \
     "$TSITSI_FILE" using 1:2 with lines lc rgb "red" lt 1 lw 3 title "Tsitsi Mean"
EOF
                
                # Run Gnuplot
                gnuplot "$GNUPLOT_SCRIPT"
                echo "Generated plot: $PLOT_FILE"
            fi
        done
    fi
done

echo "Comparison complete. Check the '$OUTPUT_DIR' directory for plots."
