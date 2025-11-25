"""
Parameter definitions for Semiconductor Bloch Equations solver.

This module contains physical constants and simulation parameters.
All energies in meV, times in fs unless otherwise specified.
"""

import numpy as np

# Physical constants
hbar = 658.5  # meV·fs 

# Laser parameters
chi0_values = [0.1]   # Coupling strength values to scan 
delta_t = 25.0         # fs (laser pulse width)
Delta_0 = 30.0         # meV (detuning - transition energy)

# Time parameters
dt = 2.0                 # fs (time step)
t_max = 1000.0          # fs (maximum time)
t_0 = -3 * delta_t       # fs (initial time)

# Energy level parameters
N = 100                  # Number of energy levels (reduced from 1000)
epsilon_max = 300.0      # meV (maximum energy level)
Delta_epsilon = epsilon_max / N  # meV (energy spacing = 3 meV)

# Dephasing parameters
T2 = 200.0                     # fs (constant dephasing time)
T2_0 = 210.0                   # fs (initial dephasing time at t=0)
gamma_scattering = 6.5e-20     # cm³/fs (scattering coefficient)

# Rydberg energy and normalization
E_R = 4.2                 # meV (Rydberg energy)
V = 1.0                   # cm³ (normalization volume)
a0 = 1.25e-6              # cm (Bohr radius)    
C0 = V * Delta_epsilon**(3/2) / (2 * np.pi**2 * a0**3 * E_R**(3/2))


# Default parameter dictionary for SBESolver
PARAMS = {
    'N': N,
    'Delta_epsilon': Delta_epsilon,
    'Delta_0': Delta_0,
    'E_R': E_R,
    'T2': T2,              
    'T2_0': T2_0,           
    'gamma': gamma_scattering, 
    'chi_0': chi0_values[0],
    'delta_t': delta_t,
    'hbar': hbar,
    't_0': t_0,
    't_max': t_max,
    'dt': dt,
    'C_0': C0,
}


