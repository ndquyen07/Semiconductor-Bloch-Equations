# Semiconductor Bloch Equations Solver

A numerical solver for the Semiconductor Bloch Equations (SBEs) describing carrier dynamics in semiconductors under optical excitation.

## Overview

The Semiconductor Bloch Equations describe the quantum dynamics of electron-hole pairs in semiconductors interacting with light. This implementation solves:

### The Equations

**Microscopic Polarization:**
```
∂p_k/∂t = -iω_k p_k - iΩ(t)(1 - n_k^e - n_k^h) - γ p_k
```

**Carrier Populations:**
```
∂n_k^e/∂t = 2 Im[Ω*(t) p_k] - γ_e n_k^e
∂n_k^h/∂t = 2 Im[Ω*(t) p_k] - γ_h n_k^h
```

Where:
- `p_k` - Microscopic polarization (coherence)
- `n_k^e, n_k^h` - Electron and hole occupation numbers
- `ω_k` - Transition energy
- `Ω(t) = μ_k E(t)` - Rabi frequency
- `γ, γ_e, γ_h` - Dephasing and relaxation rates

## Features

- ✓ Multiple optical field types (Gaussian, CW, chirped, sech)
- ✓ Pauli blocking (phase space filling)
- ✓ Dephasing and relaxation processes
- ✓ Optional Coulomb interactions (Hartree-Fock)
- ✓ Comprehensive visualization
- ✓ Emission spectrum calculation

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic Example

```python
python main.py
```

### Custom Parameters

```python
from src.sbe_solver import SBESolver
from src.fields import GaussianPulse
import numpy as np
from scipy.integrate import solve_ivp

# Define parameters
params = {
    'omega_k': 2.0,      # Transition energy (eV)
    'gamma': 0.1,        # Dephasing rate
    'gamma_e': 0.05,     # Electron relaxation
    'gamma_h': 0.05,     # Hole relaxation
    'mu_k': 1.0,         # Dipole moment
}

# Create optical field
field = GaussianPulse(amplitude=0.5, t0=20.0, sigma=5.0, omega=2.0)

# Initialize solver
solver = SBESolver(params, field)

# Solve
solution = solve_ivp(
    solver.derivatives,
    (0, 100),
    [0, 0, 0, 0],  # Initial state
    method='RK45',
    dense_output=True
)
```

## Physical Interpretation

### Key Phenomena

1. **Rabi Oscillations**: Coherent oscillations between ground and excited states
2. **Dephasing**: Loss of quantum coherence (γ term)
3. **Pauli Blocking**: Reduction of transition probability when states are occupied (1 - n_e - n_h factor)
4. **Carrier Generation**: Optical excitation creates electron-hole pairs
5. **Relaxation**: Population decay back to equilibrium

### Typical Parameter Ranges

- **Transition energy**: 1-3 eV (visible to near-IR)
- **Dephasing time**: 10-100 fs (1/γ)
- **Relaxation time**: 100-1000 fs (1/γ_e)
- **Pulse duration**: 5-50 fs (ultrafast regime)

## Numerical Methods

The solver uses:
- **ODE Integration**: Runge-Kutta 4/5 (RK45) from scipy
- **Complex Variables**: Split into real/imaginary parts for real ODE solver
- **Adaptive Stepping**: Automatic time step adjustment for accuracy

## Advanced Features

### Including Coulomb Interactions

```python
from src.sbe_solver import SBESolverWithCoulomb

# Coulomb matrix element
V_kk = 0.1  # Interaction strength

solver = SBESolverWithCoulomb(params, field, V_kk)
```

### Different Field Types

```python
from src.fields import CWField, ChirpedPulse, SechPulse

# Continuous wave
field = CWField(amplitude=0.3, omega=2.0)

# Chirped pulse
field = ChirpedPulse(amplitude=0.5, t0=20, sigma=5, omega=2.0, chirp=0.01)

# Sech pulse (ultrafast)
field = SechPulse(amplitude=0.5, t0=20, tau=3.0, omega=2.0)
```

## Output

The solver generates:
- Time evolution plots of polarization and populations
- Phase space trajectories
- External field visualization
- PNG figures saved in the working directory

## References

1. Haug, H., & Koch, S. W. (2009). *Quantum Theory of the Optical and Electronic Properties of Semiconductors*
2. Kira, M., & Koch, S. W. (2011). *Semiconductor Quantum Optics*
3. Schäfer, W., & Wegener, M. (2002). *Semiconductor Optics and Transport Phenomena*

## License

MIT License
