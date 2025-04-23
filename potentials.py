import numpy as np
import matplotlib.pyplot as plt

# Define x as a 1D field (radial direction)
x = np.linspace(-2, 2, 400)

# Define potential for three cases of mu^2
lambda_val = 0.5

# Case 1: mu^2 < 0
mu_sq_neg = -1
V_neg = mu_sq_neg * x**2 + lambda_val * x**4

# Case 2: mu^2 = 0
mu_sq_zero = 0
V_zero = mu_sq_zero * x**2 + lambda_val * x**4

# Case 3: mu^2 > 0
mu_sq_pos = 1
V_pos = mu_sq_pos * x**2 + lambda_val * x**4

# Plotting
plt.figure(figsize=(5, 4))
plt.plot(x, V_neg, label=r'$m_{HH}^2 < 0$', color='red')
plt.plot(x, V_zero, label=r'$m_{HH}^2 = 0$', color='green')
plt.plot(x, V_pos, label=r'$m_{HH}^2 > 0$', color='blue')
plt.grid()
plt.legend(fontsize=14)
plt.title(r"Higgs Potential $V(x) = m_{HH}^2 x^2 + \lambda_{1} x^4$",fontsize=14)
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$V(x)$",fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()

plt.legend()
plt.savefig('higgs_potentials.pdf', format='pdf')
