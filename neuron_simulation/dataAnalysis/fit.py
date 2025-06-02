import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf
import os

# Load data
filename = "MSDisp_statistics.txt"
base_directory = os.path.join("..", "simulations", "Ising_0,9_10,0_6,0", "5k_cluster_100Lx_0,5rad_0,39270pac_10e5time_1branch_2synapses_statistics")
file_path = os.path.join(base_directory, filename)
data = np.genfromtxt(file_path, skip_header=1)
x_data = data[:, 0]
y_data = data[:, 1]

# Define fitting functions
def power_law(x, A, B):
    return B * x ** A

def modified_function(x, A, B):
    return A*(B * x) / (1 + B*x)

def error_function(x, A, B, C):
    return A * erf(B * (x - C))

# Fit all three models
try:
    # Power law fit
    p0_power = [0, 0.1]
    params_power, _ = curve_fit(power_law, x_data, y_data, p0=p0_power)
    A_pow, B_pow = params_power
    
    # Modified function fit
    p0_mod = [1, 0.001]
    params_mod, _ = curve_fit(modified_function, x_data, y_data, p0=p0_mod)
    A_mod, B_mod = params_mod
    
    # Error function fit (with bounds to prevent unstable fits)
    p0_erf = [max(y_data), 0.001, np.median(x_data)]
    params_erf, _ = curve_fit(error_function, x_data, y_data, p0=p0_erf, 
                             bounds=([0, 0, 0], [2*max(y_data), 0.1, max(x_data)]))
    A_erf, B_erf, C_erf = params_erf
    
except RuntimeError as e:
    print(f"Curve fitting failed: {e}")

# Generate smooth x values for plotting curves
x_smooth = np.linspace(min(x_data), max(x_data), 500)

# Plot all fits
plt.figure(figsize=(12, 8))
plt.scatter(x_data, y_data, color='black', label="Experimental Data", zorder=1)

# Plot power law fit
plt.plot(x_smooth, power_law(x_smooth, A_pow, B_pow), 
         color='blue', linestyle='--', 
         label=f"Power Law: y = {B_pow:.3f}x^{A_pow:.3f}")

# Plot modified function fit
plt.plot(x_smooth, modified_function(x_smooth, A_mod, B_mod), 
         color='green', linestyle='-.', 
         label=f"Modified: y = ({A_mod:.3f}x)/({B_mod:.3f}+x)")

# Plot error function fit
plt.plot(x_smooth, error_function(x_smooth, A_erf, B_erf, C_erf), 
         color='red', linestyle=':', 
         label=f"Error Func: y = {A_erf:.3f}·erf({B_erf:.5f}(x-{C_erf:.1f}))")

plt.xlabel("Time Steps", fontsize=12)
plt.ylabel("Connected Fraction", fontsize=12)
plt.title("Comparison of Fitting Functions", fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Print fitting results
print("\nFitting Results:")
print(f"Power Law: A = {A_pow:.5f}, B = {B_pow:.5f}")
print(f"Modified Function: A = {A_mod:.5f}, B = {B_mod:.5f}")
print(f"Error Function: A = {A_erf:.5f}, B = {B_erf:.5f}, C = {C_erf:.5f}")

plt.show()