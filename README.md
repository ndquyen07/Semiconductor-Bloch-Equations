# Semiconductor Bloch Equations Solver

A numerical solver for the multi-level Semiconductor Bloch Equations (SBEs) describing carrier dynamics in semiconductors under optical excitation with Coulomb many-body interactions.

## Overview

This implementation solves the multi-level Semiconductor Bloch Equations for GaAs-like semiconductors, incorporating:
- Multi-level energy structure (N discrete levels)
- Coulomb many-body interactions via Hartree-Fock approximation
- Time-dependent dephasing due to carrier-carrier scattering
- Gaussian laser pulse excitation
- Absorption spectrum calculation via Fourier transform

### The Equations

**Carrier Population Dynamics:**
$$\frac{\partial f_n}{\partial t} = -2\text{Im}\left[\Omega_n^R(t) p_n^*(t)\right]$$

**Microscopic Polarization:**
$$\frac{\partial p_n}{\partial t} = -\frac{i}{\hbar}\left[n\Delta\varepsilon - \Delta_0 - E_n(t)\right]p_n + i\left[1 - 2f_n(t)\right]\Omega_n^R(t) - \frac{p_n(t)}{T_2(t)}$$

**Rabi Frequency with Coulomb Interactions:**
$$\Omega_n^R(t) = \chi(t) + \frac{1}{\hbar}\frac{\sqrt{E_R}}{\pi}\Delta\varepsilon \sum_{n'=1}^N g(n,n')p_{n'}(t)$$

**Coulomb Energy Shift:**
$$E_n(t) = \frac{\sqrt{E_R}}{\pi}\Delta\varepsilon \sum_{n'=1}^N g(n,n')2f_{n'}(t)$$

**Coupling Matrix:**
$$g(n,n') = \frac{1}{\sqrt{n\Delta\varepsilon}}\ln\left|\frac{\sqrt{n} + \sqrt{n'}}{\sqrt{n} - \sqrt{n'}}\right|$$

**Time-Dependent Dephasing:**
$$\frac{1}{T_2(t)} = \frac{1}{T_{2,0}} + \gamma N(t)$$

Where:
- $f_n(t)$ - Occupation probability of level $n$
- $p_n(t)$ - Microscopic polarization at level $n$ (complex)
- $\Omega_n^R(t)$ - Effective Rabi frequency (includes Coulomb field)
- $E_n(t)$ - Coulomb-induced energy shift
- $\chi(t)$ - Gaussian laser pulse
- $\Delta\varepsilon$ - Energy spacing between levels
- $\Delta_0$ - Detuning from band edge
- $E_R$ - Rydberg energy (4.2 meV for GaAs)
- $T_2(t)$ - Time-dependent dephasing time
- $\gamma$ - Scattering coefficient
- $N(t)$ - Total carrier density

## Features

- ✅ Multi-level quantum system (up to 100+ energy levels)
- ✅ Full Coulomb many-body interactions (Hartree-Fock)
- ✅ Carrier-density-dependent dephasing
- ✅ Pauli blocking (phase space filling)
- ✅ Gaussian laser pulse excitation
- ✅ Absorption spectrum calculation
- ✅ Custom 4th-order Runge-Kutta solver
- ✅ Comparison mode (with/without Coulomb interactions)
- ✅ Comprehensive visualization tools

## Installation

```bash
pip install -r requirements.txt
```

**Requirements:**
- Python 3.11+
- NumPy ≥ 1.21.0
- Matplotlib ≥ 3.4.0

## Usage

### Quick Start

Run the default simulation for GaAs:

```bash
python main.py
```

This will:
1. Solve SBEs with and without Coulomb interactions
2. Save results to `data/` directory
3. Generate plots in `result/` directory
4. Compare absorption spectra

### Custom Simulation

```python
from src.sbe_solver import SBESolver
from src.params import PARAMS

# Create solver instance
solver = SBESolver()

# Run simulation with custom parameters
solver.fit(
    chi_0=0.001,         # Laser coupling strength
    delta_t=25.0,        # Laser pulse width (fs)
    Delta_0=30.0,        # Detuning (meV)
    E_R=4.2,             # Rydberg energy (meV)
    T2_0=210.0,          # Initial dephasing time (fs)
    gamma=6.5e-20,       # Scattering coefficient (cm³/fs)
    hbar=658.5,          # ℏ in meV·fs
    N=100,               # Number of energy levels
    Delta_epsilon=3.0,   # Energy spacing (meV)
    t_0=-75.0,           # Start time (fs)
    t_max=1000.0,        # End time (fs)
    dt=2.0,              # Time step (fs)
    C_0=1.0,             # Normalization constant
    with_coulomb=True    # Include Coulomb interactions
)

# Extract results
t, f_n, p_n, population, polarization, absorption, energy = solver.extract_results()
```

### Parameter Configuration

Modify parameters in `src/params.py`:

```python
# Laser parameters
chi0 = 0.001          # Coupling strength
delta_t = 25.0        # Pulse width (fs)
Delta_0 = 30.0        # Detuning (meV)

# Energy levels
N = 100               # Number of levels
epsilon_max = 300.0   # Maximum energy (meV)

# Dephasing
T2_0 = 210.0          # Initial T2 (fs)
gamma = 6.5e-20       # Scattering (cm³/fs)

# Rydberg parameters
E_R = 4.2             # Rydberg energy (meV) for GaAs
```

## Physical Interpretation

### Key Many-Body Effects

1. **Coulomb Enhancement**: Many-body Coulomb field enhances the optical response
2. **Excitation-Induced Dephasing (EID)**: Dephasing rate increases with carrier density
3. **Band Gap Renormalization**: Coulomb interactions shift transition energies
4. **Pauli Blocking**: Occupied states reduce transition probabilities
5. **Absorption Saturation**: Nonlinear absorption at high intensities

### Parameter Ranges (GaAs)

| Parameter | Symbol | Typical Value | Unit |
|-----------|--------|---------------|------|
| Rydberg energy | $E_R$ | 4.2 | meV |
| Dephasing time | $T_{2,0}$ | 50-500 | fs |
| Pulse width | $\delta_t$ | 10-100 | fs |
| Scattering coeff. | $\gamma$ | 10⁻²⁰-10⁻¹⁹ | cm³/fs |
| Energy spacing | $\Delta\varepsilon$ | 1-10 | meV |

## Numerical Methods

### ODE Integration
- **Algorithm**: 4th-order Runge-Kutta (RK4)
- **Implementation**: Custom solver in `src/rk4_solver.py`
- **Time step**: Fixed step size (typically 2 fs)
- **Stability**: Verified for complex-valued equations

### Vectorization
- Precomputed coupling matrix $g(n,n')$ for efficiency
- Vectorized operations over all energy levels
- Matrix multiplication for Coulomb sums

### Fourier Transform
- Direct numerical integration method
- Computes absorption spectrum: $\alpha(\omega) \propto \text{Im}\left[\frac{P(\omega)}{E(\omega)}\right]$
- Energy range: -100 to 100 meV

## Project Structure

```
SemiconductorBlochEquations/
├── main.py                    # Main simulation script
├── src/
│   ├── sbe_solver.py          # Core SBE solver class
│   ├── rk4_solver.py          # RK4 integration routine
│   ├── params.py              # Physical parameters
│   ├── utils.py               # Data I/O utilities
│   └── visualization.py       # Plotting functions
├── experiment/
│   ├── experiment.ipynb       # Interactive notebooks
│   └── experiment_3d.ipynb    # 3D visualization
├── data/                      # Saved numerical results
├── result/                    # Generated plots
└── requirements.txt           # Python dependencies
```

## Output

The solver generates:

### Data Files (in `data/`)
- Time evolution of populations and polarizations
- Absorption spectra (with/without Coulomb)
- Complete state vectors

### Plots (in `result/`)
- Time-dependent population $N(t)$
- Microscopic polarization $|P(t)|$
- Laser pulse envelope $\chi(t)$
- Absorption spectrum comparison
- Level-resolved dynamics

## Examples

### Compare Coulomb Effects

```python
from src.sbe_solver import SBESolver
from src.params import PARAMS

# Solve with Coulomb interactions
solver_coulomb = SBESolver()
solver_coulomb.fit(**PARAMS, with_coulomb=True)

# Solve without Coulomb interactions
solver_free = SBESolver()
solver_free.fit(**PARAMS, with_coulomb=False)

# Compare absorption spectra
from src.visualization import plot_absorption_comparison
plot_absorption_comparison(solver_coulomb, solver_free)
```

### Access Level-Resolved Data

```python
# Get occupation of specific level (e.g., n=10)
f_10 = solver.f_n[:, 9]  # Zero-indexed

# Get polarization at all levels
p_all = solver.p_n  # Shape: (time_points, N)

# Total population density
N_t = solver.population
```

## Physics Background

This implementation is based on the multi-level SBE formalism for modeling ultrafast optical processes in semiconductors. The key physics includes:

- **Quasi-2D confinement**: Energy levels discretized in quantum well structures
- **Coulomb screening**: Static Hartree-Fock approximation
- **Dephasing mechanisms**: Carrier-carrier scattering (T₂ depends on density)
- **Sum rules**: Particle number conservation maintained numerically

## References

1. **Haug, H., & Koch, S. W.** (2009). *Quantum Theory of the Optical and Electronic Properties of Semiconductors* (5th ed.). World Scientific.
2. **Kira, M., & Koch, S. W.** (2011). *Semiconductor Quantum Optics*. Cambridge University Press.
3. **Schäfer, W., & Wegener, M.** (2002). *Semiconductor Optics and Transport Phenomena*. Springer.
4. **Lindberg, M., & Koch, S. W.** (1988). Effective Bloch equations for semiconductors. *Physical Review B*, 38(5), 3342.

## License

MIT License

## Author

Nguyen Dinh Quyen  
Vietnam National University - Ho Chi Minh City University of Science
