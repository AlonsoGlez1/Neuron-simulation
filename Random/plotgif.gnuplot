# Set the output format and file name
set terminal gif animate size 800,800 delay 100
set output 'prueba1.gif'

# Set the plot range and ratio
set size square
set xrange [0:500]
set yrange [0:500]

# Set up the plot style
set style data lines
set style line 1 lc rgb 'black' pt 10 ps 1.5
set style line 2 lc rgb 'black' lt 1 lw 2

# Plot circles once
plot 'particles' using 1:2:3 with circles lc rgb "blue" fill solid 0.2 notitle

# Loop through the lines
do for [i=1:1000] {
    # Replot circles to keep them in the plot
    plot 'particles' using 1:2:3 with circles lc rgb "blue" fill solid 0.2 notitle, 'connections' index 0 every ::0::i with lines ls 2 notitle

    # Pause to create frames for the GIF
    pause 0.05
}
