from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd

file_path = 'scan.dat.gz'

data = pd.read_csv(file_path, sep=r'\s+', low_memory= False)

alpha = 1/137 #fine structure constant,  = e**2/(4 * pi * epsilon_0 * h_bar * c)
nu = 246 #VEV

def f_a(x):
    return -5 + 12 * np.log(x)

def f_b(x):
    return 3 - 4 * np.log(x)

def f_c(x, y):
    mask = np.isclose(x, y, rtol=1e-10)
    result = np.zeros_like(x)
    result[~mask] = ((x[~mask] + y[~mask]) / 2) - ((x[~mask] * y[~mask]) / (x[~mask] - y[~mask])) * np.log(x[~mask] / y[~mask])
    return result

data['MDP'] = data['DMP'] + data['MD1']
data['MD2'] = data['DM3'] + data['DMP'] + data['MD1']
data['DM2'] = data['DM3'] + data['DMP']
    
data['x1'] = data['MD1']/data['MDP']
data['x2'] = data['MD2']/data['MDP']

data["S"] = ((1/(72*np.pi)) * (1/((data['x2']**2 - data['x1']**2)**3))
            * (((data['x2']**6) * f_a(data['x2'])) - ((data['x1']**6)
            * f_a(data['x1'])) + (9 * (data['x2']**2) * (data['x1']**2 )
            * (((data['x2']**2) * f_b(data['x2'])) - ((data['x1']**2)
                                                                        * f_b(data['x1']))))))
data["T"] = ((1/(32*(np.pi**2)*alpha*(nu**2)))
            * (f_c(data["MDP"]**2,data["MD2"]**2) + f_c(data["MDP"]**2,data["MD1"]**2) - f_c(data["MD2"]**2,data["MD1"]**2)))

version = 'all'

S_central, S_error = 0.06, 0.09
T_central, T_error = 0.1, 0.07



"""
if version == 'less':
    cutT = (data['T'] < (T_central + T_error))
    cutS = (data['S'] < (S_central + S_error))

if version == 'more':
    cutT = (data['T'] > (T_central - T_error))
    cutS = (data['S'] > (S_central - S_error))

if version == 'all':
    cutT = (data['T'] > (T_central - T_error)) & (data['T'] < (T_central + T_error))
    cutS = (data['S'] > (S_central - S_error)) & (data['S'] < (S_central + S_error))

cuts = cutT & cutS
data = data[cuts]
    
cuts = cutT & cutS

data = data[cuts]
"""
print(np.shape(data))
#plt.scatter(data['S'], data['T'], s = 0.1)

#plt.close('all')

plt.axvline(x = S_central - S_error, color = "black", lw = 0.1)
plt.axvline(x = S_central + S_error, color = "black", lw = 0.1)
plt.axhline(y = T_central + T_error, color = "black", lw = 0.1)
plt.axhline(y = T_central - T_error, color = "black", lw = 0.1)
 
plt.scatter(data['S'], data['T'], s=0.5)

from matplotlib.patches import Ellipse

# Experimental values from PDG (central values and standard deviations)
S_central = 0.05
T_central = 0.09
sigma_S = 0.10
sigma_T = 0.13

# Define chi-squared thresholds for 68% and 95% confidence levels
chi2_68 = 2.30  # 68% confidence
chi2_95 = 6.18  # 95% confidence

# Create mesh grid for plotting
S_vals = np.linspace(S_central - 3*sigma_S, S_central + 3*sigma_S, 200)
T_vals = np.linspace(T_central - 3*sigma_T, T_central + 3*sigma_T, 200)
S_grid, T_grid = np.meshgrid(S_vals, T_vals)

# Compute the chi-squared function (Gaussian form for the ellipse)
def chi_squared(S, T):
    return ((S - S_central)**2 / sigma_S**2) + ((T - T_central)**2 / sigma_T**2)

# Calculate chi-squared values
chi2_vals = chi_squared(S_grid, T_grid)

# Plot the contours for 68% and 95% confidence regions
fig, ax = plt.subplots(figsize=(8, 6))
CS = ax.contour(S_grid, T_grid, chi2_vals, levels=[chi2_68, chi2_95], colors=['blue', 'red'], linewidths=2)

# Add labels and plot configuration
ax.set_xlabel('S')
ax.set_ylabel('T')
ax.set_title('S-T Contour Plot (Experimental Confidence Regions)')

# Add legends for the contour levels
ax.clabel(CS, inline=True, fontsize=12, fmt='%.2f')


plt.xlim(-0.2, 0.3)
plt.ylim(-0.1, 0.3)
plt.xlabel('S')
plt.ylabel('T')
plt.legend()
plt.show()