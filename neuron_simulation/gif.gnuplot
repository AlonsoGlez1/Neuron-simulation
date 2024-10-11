# Set the output format and file name
set terminal gif animate size 1200,1200 delay 100
set output '500.gif'

# Set the plot range and ratio
set size square
set xrange [0:50]
set yrange [0:50]

# Set margins to create space for the label
set lmargin at screen 0.1
set rmargin at screen 0.9
set tmargin at screen 0.9
set bmargin at screen 0.1

# Set up the plot style
set style data lines
set style line 1 lc rgb "black" lw 2              # Line style 1: 
set style line 2 lc rgb "dark-red" lw 2           # Red connections
set style line 3 lc rgb "dark-spring-green" lw 2  # Green connections
set style line 4 lc rgb "orchid4" lw 2            # Purple connections

# Define array of colors for circles
colors = "red green purple"

# Plot all neurons by themselves and label for the current time step
set label 1 sprintf("Time Step: %d", 0) at screen 0.1, 0.95 font "Times-Roman,24" tc rgb "black"
plot 'neurons_dat' using 1:2:3 with circles lc rgb "blue" fill solid 0.2 notitle
    
# Loop through the connections and plot them 
do for [i=1:50] {
    # Update label for the current time step
    set label 1 sprintf("Time Step: %d", i) at screen 0.1, 0.95 font "Times-Roman,24" tc rgb "black"
   
    # Replot neurons to keep them in the plot and plotting the connecting neurons with colors
    plot 'neurons_dat' using 1:2:3 with circles lc rgb "blue" fill solid 0.2 notitle, \
         for [j=0:2] 'connections_dat' index j using 1:2:3 every ::1::i with circles lc rgb word(colors, j+1) fill solid 0.5 notitle, \
         for [j=0:2] 'connections_dat' index j every ::1::i with lines ls (j+2) notitle, \
         'connections_dat' index 0 using 1:2:3 every ::1::1 with circles lc rgb "black" fill solid 1 notitle

    # Pause to create frames for the GIF
    pause 0.05
}

