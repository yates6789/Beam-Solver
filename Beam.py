# Beam.py
# James Yates - 24/05/2020
# This script implements the static beam equation to a cantilevered beam

# QUICK TEST

# Nomenclature ---------------------------------------------------------------------------------------------------------
# h - Beam Height - [m]
# b - Beam Width - [m]
# L - Beam Length - [m]
# x - Local Position - [m]
# I - Second Moment of Area - [mm4]
# E - Young's Modulus - [Pa]
# q - Distributed Load - [N/m]
# V_p - Magnitude of Point Loads - [N]
# L_p - Position of Point Loads - [m]
# V - Shear Force - [N]
# M - Bending Moment - [N*m]
# Y - Deflection - [m]

# Importing Distributions ----------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# Allowable Error ------------------------------------------------------------------------------------------------------
Error = 1e-10

# Beam Dimensions ------------------------------------------------------------------------------------------------------
h = 0.01
b = 0.01
L = 1

# Moment of Area -------------------------------------------------------------------------------------------------------
x = np.linspace(0, L, 10)
I = ((b * h ** 3) / 12) * (np.linspace(1, 1, x.size))

# Loading --------------------------------------------------------------------------------------------------------------
E = 180e9
q = -np.linspace(0, 0, x.size)
V_p = np.array([5, 2.5])
L_p = np.array([0.33, 0.66])

# Shear Force - Distributed Load - Finite Difference Scheme ------------------------------------------------------------
V = np.zeros(x.size)
ea_V = 1
if sum(q == 0) != x.size:
    while ea_V > Error:
        V_1 = V.sum()
        V[0:-1] = V[1:] - q[0:-1] * (x[1] - x[0])
        ea_V = (V.sum() - V_1) / (V.sum())
        ea_V = np.abs(ea_V)

# Shear Force - Point Loads --------------------------------------------------------------------------------------------
for i in range(0, V_p.size):
    V[x > L_p[i]] = V[x > L_p[i]] - V_p[i]
V = V - V[-1]

# Bending Moment - Finite Difference Scheme ----------------------------------------------------------------------------
M = np.zeros(x.size)
ea_M = 1
while ea_M > Error:
    M_1 = M.sum()
    M[0:-1] = M[1:] + V[0:-1] * (x[1] - x[0])
    ea_M = (M.sum() - M_1) / (M.sum())
    ea_M = np.abs(ea_M)

# Deflection - Finite Difference Scheme --------------------------------------------------------------------------------
Y = np.zeros(x.size)
ea_Y = 1
while ea_Y > Error:
    Y_1 = Y.sum()
    Y[2:] = 2 * Y[1:-1] - Y[0:-2] - (M[2:] / (E * I[2:])) * (x[1] - x[0]) ** 2
    ea_Y = (Y.sum() - Y_1) / (Y.sum())
    ea_Y = np.abs(ea_Y)

# Plots ----------------------------------------------------------------------------------------------------------------
fig = plt.figure()
# Distributed Load
ax_q = fig.add_subplot(2, 2, 1)
ax_q.plot(x, q)
ax_q.set(title='Distributed Load', xlabel='x - [m]', ylabel='q(x) - [N/m]')
plt.grid(True)
# Shear Force
ax_V = fig.add_subplot(2, 2, 2)
ax_V.plot(x, V)
ax_V.set(title='Shear Force', xlabel='x - [m]', ylabel='V(x) - [N]')
plt.grid(True)
# Bending Moment
ax_M = fig.add_subplot(2, 2, 3)
ax_M.plot(x, M)
ax_M.set(title='Bending Moment', xlabel='x - [m]', ylabel='M(x) - [N*m]')
plt.grid(True)
# Deflection
ax_M = fig.add_subplot(2, 2, 4)
ax_M.plot(x, Y)
ax_M.set(title='Deflection', xlabel='x - [m]', ylabel='Y(x) - [m]')
plt.grid(True)
# Printing Plot
fig.subplots_adjust(wspace=0.3, hspace=0.5)
plt.show()
