import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# Experimental values from PDG (central values and standard deviations)
S_central = -0.01
T_central = 0.04
sigma_S = 0.07
sigma_T = 0.06

# Assuming a correlation coefficient between S and T (from PDG review)
# If the correlation coefficient is given explicitly in the paper, replace the value.
Corr_ST = 0.92  # Example correlation, adjust based on PDG review

# Covariance matrix
Cov_ST = Corr_ST * sigma_S * sigma_T
cov_matrix = np.array([[sigma_S**2, Cov_ST], [Cov_ST, sigma_T**2]])

# Define chi-squared thresholds for 68% and 95% confidence levels
chi2_68 = 2.30  # 68% confidence
chi2_95 = 6.18  # 95% confidence

# Calculate the angle of the ellipse (in radians)
angle = 0.5 * np.arctan2(2 * Cov_ST, sigma_S**2 - sigma_T**2)  # in radians
angle_deg = np.degrees(angle)  # convert to degrees for plotting

# Create mesh grid for plotting
S_vals = np.linspace(S_central - 3*sigma_S, S_central + 3*sigma_S, 200)
T_vals = np.linspace(T_central - 3*sigma_T, T_central + 3*sigma_T, 200)
S_grid, T_grid = np.meshgrid(S_vals, T_vals)

# Compute the chi-squared function (Gaussian form for the ellipse)
def chi_squared(S, T):
    diff = np.array([S - S_central, T - T_central])
    return np.dot(diff.T, np.linalg.solve(cov_matrix, diff))

# Calculate chi-squared values
chi2_vals = np.vectorize(chi_squared)(S_grid, T_grid)

# Plotting only the ellipses at chi-squared values for 68% and 95% confidence regions
fig, ax = plt.subplots(figsize=(8, 6))

# Create Ellipse patches for 68% and 95% contours
ellipse_68 = Ellipse((S_central, T_central), width=2 * sigma_S * np.sqrt(chi2_68), height=2 * sigma_T * np.sqrt(chi2_68), angle=angle_deg, edgecolor='blue', facecolor='none', linewidth=2)
ellipse_95 = Ellipse((S_central, T_central), width=2 * sigma_S * np.sqrt(chi2_95), height=2 * sigma_T * np.sqrt(chi2_95), angle=angle_deg, edgecolor='red', facecolor='none', linewidth=2)

# Add the ellipses to the plot
ax.add_patch(ellipse_68)
ax.add_patch(ellipse_95)

# Add labels and plot configuration
ax.set_xlabel('S')
ax.set_ylabel('T')
ax.set_title('S-T Contour Plot (Tilted Ellipses for Confidence Regions)')
ax.set_xlim(-1, 1)
ax.set_ylim(-0.8, 0.8)

# Display the plot
plt.grid(True)
plt.show()

# Print the angle of the ellipse for reference
print(f"The angle of the ellipse is {angle_deg:.2f} degrees.")
