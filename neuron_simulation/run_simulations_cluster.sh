#!/bin/bash

# Fixed arguments
ARG1=125
ARG2=100000
ARG3="neurons_dat_5k_125Lx_cluster"
BASE_FILENAME_CONNECTIONS="connections_dat"
BASE_FILENAME_RESULTS="results_dat"
OUTPUT_DIR1="simulations"
OUTPUT_DIR2="5k_cluster_125Lx_0,5rad_0,25133pac_10e5time_1branch_2synapses"


# Loop from 1 to 50
for i in {1..50}; do
    # Create the full filename with the current value of i
    FILENAME1="${OUTPUT_DIR1}/${OUTPUT_DIR2}/${BASE_FILENAME_CONNECTIONS}_${i}"
    FILENAME2="${OUTPUT_DIR1}/${OUTPUT_DIR2}/${BASE_FILENAME_RESULTS}_${i}"


    # Display the current command
    echo "Executing: ./runNeurons $ARG1 $ARG2 $ARG3 $FILENAME1 $FILENAME2"

    # Run the program
    ./runNeurons "$ARG1" "$ARG2" "$ARG3" "$FILENAME1" "$FILENAME2"

    # Check the exit status
    if [ $? -ne 0 ]; then
        echo "Error: Program failed on iteration $i with filename $FILENAME1" 
        exit 1
    fi
done

echo "All executions completed successfully!"
