#!/bin/bash

run_simulation() {
    local sim_type=$1
    local box_size=$2
    local pac_value=$3

    local positions_path="positions/neurons_dat_5k_${box_size}Lx_$sim_type"
    local output_dir="5k_${sim_type}_${box_size}Lx_0,5rad_${pac_value}pac_10e5time_1branch_2synapses"

    echo "Running ${sim_type^} Simulation (${box_size}Lx)..."

    mkdir -p "simulations/${output_dir}" 2>/dev/null

    for i in {1..50}; do
        local conn_file="simulations/${output_dir}/connections_dat_${i}"
        local res_file="simulations/${output_dir}/results_dat_${i}"

        echo "Executing: ./runNeurons $box_size 100000 $positions_path $conn_file $res_file"
        
        if ! ./runNeurons "$box_size" 100000 "$positions_path" "$conn_file" "$res_file"; then
            echo "Error: Program failed on iteration $i with filename $conn_file" 
            return 1
        fi
    done

    echo "${sim_type^} simulation (${box_size}Lx) completed successfully!"
    return 0
}

# Parameter pairs: box_size and pac_value
declare -A params=(
    [100]="0,39270"
    [125]="0,25133"
    [150]="0,17453"
    [175]="0,12823"
    [200]="0,09817"
)

# Run all simulations for each parameter pair and each type
for box_size in "${!params[@]}"; do
    pac_value="${params[$box_size]}"
    for sim_type in random lattice hexagonal cluster; do
        if ! run_simulation "$sim_type" "$box_size" "$pac_value"; then
            echo "${sim_type^} simulation failed for ${box_size}Lx. Exiting."
            exit 1
        fi
    done
done

echo "All simulations completed successfully!"