import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from scipy.spatial.distance import mahalanobis

# Given data: means, standard deviations, and correlation coefficient
mu_1, mu_2 = 0.06, 0.1  # Means for S and T
sigma_1, sigma_2 = 0.09, 0.07  # Standard deviations for S and T
rho = 0.91  # Correlation coefficient

# Covariance matrix
cov_matrix = np.array([
    [sigma_1**2, rho * sigma_1 * sigma_2],
    [rho * sigma_1 * sigma_2, sigma_2**2]
])

# Inverse of covariance matrix for Mahalanobis distance calculation
cov_inv = np.linalg.inv(cov_matrix)

# Generate random points within a range
n_points = 1000
np.random.seed(42)
random_points = np.random.multivariate_normal([mu_1, mu_2], cov_matrix, size=n_points)
print(random_points[1])

# Calculate Mahalanobis distances and filter points within 2 standard deviations
distances = np.array([mahalanobis(point, [mu_1, mu_2], cov_inv) for point in random_points])
filtered_points = random_points[distances <= np.sqrt(chi2.ppf(0.95, df=2))]
unfiltered_points = random_points[distances > np.sqrt(chi2.ppf(0.95, df=2))]

# Plot the ellipses without data points
fig, ax = plt.subplots(figsize=(8, 6))

# Confidence ellipse plotting function
def confidence_ellipse(mean, cov, ax, n_std=1, **kwargs):
    eigvals, eigvecs = np.linalg.eigh(cov)
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
    scale_factor = np.sqrt(chi2.ppf({1: 0.68, 2: 0.95}[n_std], df=2))
    width, height = 2 * np.sqrt(eigvals) * scale_factor
    ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, **kwargs)
    ax.add_patch(ellipse)

# Plot ellipses for 1st and 2nd standard deviation
for ns in [1, 2]:  # Only plot 1st and 2nd standard deviation ellipses
    if ns == 1:
        color = 'red'
    else:
        color = 'blue'
    confidence_ellipse([mu_1, mu_2], cov_matrix, ax, n_std=ns, edgecolor=color, fill=False)

# Plot filtered points (within 2 standard deviations)
ax.scatter(filtered_points[:, 0], filtered_points[:, 1], color='green', alpha=0.5, label='Filtered Points (within 2 std)')

# Plot unfiltered points (outside 2 standard deviations)
ax.scatter(unfiltered_points[:, 0], unfiltered_points[:, 1], color='orange', alpha=0.5, label='Unfiltered Points (outside 2 std)')

# Set axis limits to ensure ellipses are visible
ax.set_xlim(mu_1 - 3 * sigma_1, mu_1 + 3 * sigma_1)
ax.set_ylim(mu_2 - 3 * sigma_2, mu_2 + 3 * sigma_2)
ax.set_aspect('equal')  # Keep the aspect ratio correct

ax.set_xlabel("S")
ax.set_ylabel("T")
plt.title("Confidence Ellipses for S and T (1st and 2nd Std Deviations) with Filtered and Unfiltered Points")
plt.grid(True)
plt.legend()
plt.show()
