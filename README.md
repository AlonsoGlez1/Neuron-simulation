# Neuron Simulation in a 2D Environment

This repository contains a C program that simulates neurons in a 2D Petri dish-like environment. The program models particle behavior, probabilistic connections based on distance, and other statistical properties, providing insights into simulated neural networks.

## Features

- **2D Neuron Placement:** Neurons are arranged in a 2D space.
- **Connection Probability:** Connections between neurons are based on probabilistic functions of distance and interdistance.
- **Adjustable Parameters:** Modify key simulation parameters like decay factors, scale factors, and packing fraction.
- **Efficiency Enhancements:** Optimized using proximity-based calculations for performance improvements.
- **Data Output:** Generates connection data and simulation results for visualization and analysis.

## Dependencies

This program requires the following libraries:

- Standard C libraries: `stdio.h`, `stdlib.h`, `string.h`, `time.h`, `math.h`, `stdbool.h`

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/neuron-simulation.git
   cd neuron-simulation
   ```

2. Compile the program:

   ```bash
   gcc neurons.c -o neuron_simulation -lm
   ```

3. Run the executable:

   ```bash
   ./neuron_simulation
   ```

## Configuration

Adjust the following macros in the `neurons.c` file to customize your simulation:

- `PI`: Value of pi (default: `3.141592654`)
- `DECAY_FACTOR`: Factor for distance decay in connection probability
- `INT_DIST_FACTOR`: Factor for interdistance decay
- `SCALE_FACTOR`: Scaling factor for the probability distribution
- `PACKING_FRACTION`: Fraction of area occupied by neurons

## Outputs

The program generates data files containing:

1. **Neuron Positions:** Coordinates of neurons in the 2D environment.
2. **Connection Data:** Probabilistic connections between neurons based on the simulation.

These outputs can be visualized using tools like GNUPlot or Python for deeper analysis.

## Usage Example

Here's how you can visualize the connections using GNUPlot:

```gnuplot
set term png
set output 'connections.png'
plot 'neuron_positions.dat' using 1:2 with points title 'Neurons', \
     'connections.dat' using 1:2:3:4 with vectors title 'Connections'
```

## Contribution

Contributions are welcome! Feel free to open issues or submit pull requests to enhance the simulation or improve the code structure.

## License

This project is licensed under the GPL-3.0 License. See the `LICENSE` file for details.

## Acknowledgments

- Inspiration from statistical physics and neural network modeling.
- Developed as part of research into neuron behavior in constrained environments.

---

### Author

Alonso\
Physics Graduate | Aspiring Programmer & Data Analyst

