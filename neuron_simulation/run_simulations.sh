#!/bin/bash

# Define a function for each simulation type

run_random_simulation() {
    echo "Running Random Simulation..."

    # Fixed arguments
    ARG1=175    # Box Size
    ARG2=100000 # Time steps
    ARG3="positions/neurons_dat_5k_175Lx_random"              # Positions data 
    BASE_FILENAME_CONNECTIONS="connections_dat"     # Connections fine name
    BASE_FILENAME_RESULTS="results_dat"             # Multi branch results file name
    OUTPUT_DIR1="simulations"
    OUTPUT_DIR2="5k_random_175Lx_0,5rad_0,12823pac_10e5time_1branch_2synapses"


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
}

run_lattice_simulation() {
    echo "Running Lattice Simulation..."

    # Fixed arguments
    ARG1=175    # Box size
    ARG2=100000 #Time steps
    ARG3="positions/neurons_dat_5k_175Lx_lattice"             # Positions data 
    BASE_FILENAME_CONNECTIONS="connections_dat"     # Connections fine name
    BASE_FILENAME_RESULTS="results_dat"             # Multi branch results file name
    OUTPUT_DIR1="simulations"
    OUTPUT_DIR2="5k_lattice_175Lx_0,5rad_0,12823pac_10e5time_1branch_2synapses"


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
}

run_hexagonal_simulation() {
    echo "Running Hexagonal Simulation..."

    # Fixed arguments
    ARG1=175    # Box Size
    ARG2=100000 # Time Steps
    ARG3="positions/neurons_dat_5k_175Lx_hexagonal"   # Positions data 
    BASE_FILENAME_CONNECTIONS="connections_dat"     # Connections fine name
    BASE_FILENAME_RESULTS="results_dat"             # Multi branch results file name
    OUTPUT_DIR1="simulations"
    OUTPUT_DIR2="5k_hexagonal_175Lx_0,5rad_0,12823pac_10e5time_1branch_2synapses"


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

}

run_cluster_simulation() {
    echo "Running Cluster Simulation..."

    # Fixed arguments
    ARG1=175    # Box Size
    ARG2=100000 # Time Steps
    ARG3="positions/neurons_dat_5k_175Lx_cluster"   # Positions data 
    BASE_FILENAME_CONNECTIONS="connections_dat"     # Connections fine name
    BASE_FILENAME_RESULTS="results_dat"             # Multi branch results file name
    OUTPUT_DIR1="simulations"
    OUTPUT_DIR2="5k_cluster_175Lx_0,5rad_0,12823pac_10e5time_1branch_2synapses"


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

}

# Execute each function in sequence
run_random_simulation
if [ $? -ne 0 ]; then
    echo "Random simulation failed. Exiting."
    exit 1
fi


run_lattice_simulation
if [ $? -ne 0 ]; then
    echo "Lattice simulation failed. Exiting."
    exit 1
fi

run_hexagonal_simulation
if [ $? -ne 0 ]; then
    echo "Hexagonal simulation failed. Exiting."
    exit 1
fi

run_cluster_simulation
if [ $? -ne 0 ]; then
    echo "Cluster simulation failed. Exiting."
    exit 1
fi

echo "All simulations completed successfully!"
