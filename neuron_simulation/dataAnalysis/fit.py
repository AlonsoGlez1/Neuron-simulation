import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data while skipping the header row
filename = "GyrationRadius_statistics.txt"  # Replace with your actual file
data = np.genfromtxt(filename, skip_header=1)  # Skips the first row (header)

# Extract x and y columns
x_data = data[:, 0]  # First column
y_data = data[:, 1]  # Second column

# Define the power-law function
def power_law(x, A, B):
    return B* x ** A

# Initial guess for A
initial_guess = [1, 1]  # You can change this based on expected behavior


# Perform curve fitting
params, covariance = curve_fit(power_law, x_data, y_data, p0=initial_guess)

# Extract optimized parameter
A_opt, B_opt = params
print(f"Optimized A: {A_opt}, Optimized B: {B_opt}")


# Plot the results
plt.scatter(x_data, y_data, color='purple', label="Experimental Data")  # Original data
plt.plot(x_data, power_law(x_data, A_opt, B_opt), color='green', label=f"Fit: y = x^{A_opt:.3f}")  # Fitted curve
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Power-Law Fit: y = x^A")
plt.show()

## Define the function to fit
#def modified_function(x, A, B):
#    return (A * x) / (B + x)
#
## Initial guess for A and B
#initial_guess = [1, 2500]
#
## Perform curve fitting
#params, covariance = curve_fit(modified_function, x_data, y_data, p0=initial_guess)
#
## Extract optimized parameters
#A_opt, B_opt = params
#print(f"Optimized A: {A_opt}, Optimized B: {B_opt}")
#
## Plot the results
#plt.scatter(x_data, y_data, color='purple', label="Original Data")  # Original data
#plt.plot(x_data, modified_function(x_data, A_opt, B_opt), color='green', label="Fitted Function")  # Fitted function
#plt.legend()
#plt.xlabel("x")
#plt.ylabel("y")
#plt.title("Curve Fitting")
#plt.show()
