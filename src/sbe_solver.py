"""
Core solver for Semiconductor Bloch Equations.
"""

import numpy as np
from .rk4_solver import rk4_solve

class SBESolver:
    """
    Solver for the Semiconductor Bloch Equations.
    
    The SBEs describe the dynamics of microscopic polarization and
    carrier populations in a semiconductor under optical excitation.
    """
    
    def __init__(self):
        self.chi_0 = None
        self.delta_t = None
        self.Delta_0 = None
        self.E_R = None
        self.T2 = None
        self.T2_0 = None                     # Initial dephasing time
        self.gamma = None                    # Scattering coefficient
        self.hbar = None
        self.N = None
        self.Delta_epsilon = None
        self.t_0 = None
        self.t_max = None
        self.dt = None
        self.t = None
        self.C0 = None
        self.with_coulomb = True

        self.f_n = None                      # Electron occupation
        self.p_n = None                      # Microscopic polarization
        self.population = None               # Total population
        self.polarization = None             # Total polarization
        self.absorption_spectrum = None      # Absorption spectrum
        self.spectrum_energy = None          # Energy axis for absorption spectrum
        self.T2_t = None                     # Time-dependent dephasing

        
    def precompute_g_matrix(self):
        """
        Precompute the coupling matrix g(n, n1) for all level pairs.
        """
        n_array = np.arange(1, self.N + 1)
        n_grid, n1_grid = np.meshgrid(n_array, n_array, indexing='ij')
        
        sqrt_n = np.sqrt(n_grid)
        sqrt_n1 = np.sqrt(n1_grid)
        
        numerator = sqrt_n + sqrt_n1
        denominator = sqrt_n - sqrt_n1
        
        with np.errstate(divide='ignore', invalid='ignore'):
            g_matrix = (1.0 / np.sqrt(n_grid * self.Delta_epsilon)) * np.log(np.abs(numerator / denominator))
        
        # Set diagonal to zero (n=n' gives 0/0)
        np.fill_diagonal(g_matrix, 0.0)
        
        return g_matrix
    

    def coulomb_energy(self, f_levels):
        """
        Compute energy E_n
        """
        prefactor = (np.sqrt(self.E_R) / np.pi) * self.Delta_epsilon
        E_n_array = prefactor * (self.g_matrix @ (2.0 * f_levels))
        
        return E_n_array
    

    def laser_pulse(self, t):
        """Gaussian laser pulse.
        """
        return 0.5 * (self.hbar * np.sqrt(np.pi) / self.delta_t) * self.chi_0 * np.exp(-t**2 / self.delta_t**2)
    

    def compute_T2(self, population):
        """
        Compute time-dependent dephasing time T₂(t).
        """
        T2_inv = (1.0 / self.T2_0) + self.gamma * population
        return 1.0 / T2_inv
    

    def rabi_frequency(self, t, p_levels):
        """
        Compute Rabi frequency Ω_n^R(t) for all levels using vectorized operations.
        Returns array of Ω_n^R values for n=1..N
        
        """
        prefactor = (np.sqrt(self.E_R) / np.pi) * self.Delta_epsilon
        
        coulomb_term = prefactor * (self.g_matrix @ p_levels)
        
        # Laser pulse is same for all levels
        pulse = self.laser_pulse(t)
        
        if self.with_coulomb:
            Omega_R_array = (1.0 / self.hbar) * (pulse + coulomb_term)
        else:
            Omega_R_array = (1.0 / self.hbar) * pulse
        
        return Omega_R_array
    

    def derivatives(self, t, y):
        """
        Compute the derivatives for the SBE system with N energy levels.
        """


        # Extract state variables
        f_n = y[:self.N]
        p_n = y[self.N:2*self.N] + 1j * y[2*self.N:3*self.N]


        # Compute current population and time-dependent T2
        # population = self.C0 * np.sum(np.sqrt(np.arange(1, self.N+1)) * f_n)
        # T2_current = self.compute_T2(population)

        if self.with_coulomb:
            E_n_array = self.coulomb_energy(f_n)
        else:
            E_n_array = np.zeros_like(f_n)

            
        Omega_R_array = self.rabi_frequency(t, p_n)

        # Compute derivatives
        dfdt = -2.0 * np.imag(Omega_R_array * np.conj(p_n))
        
        n_array = np.arange(1, self.N + 1)
        detuning = n_array * self.Delta_epsilon - self.Delta_0 - E_n_array
        

        dpdt = ((-1j / self.hbar) * detuning * p_n + 
                1j * (1.0 - 2.0 * f_n) * Omega_R_array - 
                p_n / self.T2)
        
        return np.concatenate([dfdt, dpdt.real, dpdt.imag])


    def fit(self, chi_0, delta_t, Delta_0, E_R, T2, T2_0, gamma, hbar, N, Delta_epsilon, t_0, t_max, dt, C_0, with_coulomb=True):
        self.chi_0 = chi_0
        self.delta_t = delta_t
        self.Delta_0 = Delta_0
        self.E_R = E_R
        self.T2 = T2
        self.T2_0 = T2_0    
        self.gamma = gamma  
        self.hbar = hbar
        self.N = N
        self.Delta_epsilon = Delta_epsilon
        self.t_0 = t_0
        self.t_max = t_max
        self.dt = dt
        self.C0 = C_0
        self.with_coulomb = with_coulomb

        # Precompute coupling matrix
        self.g_matrix = self.precompute_g_matrix()

        
        y0 = np.zeros(3 * self.N)

        # Time grid
        t_points = int((self.t_max - self.t_0) / self.dt) + 1
        t_span = (self.t_0, self.t_max)
        t_eval = np.linspace(self.t_0, self.t_max, t_points)
        

        # Solve the SBEs using RK4
        solution = rk4_solve(
            self.derivatives,
            t_span,
            y0,
            t_eval=t_eval
        )
        

        # Extract results
        self.t = solution['t']
        self.f_n = solution['y'][:self.N, :]
        self.p_n = solution['y'][self.N:2*self.N, :] + 1j * solution['y'][2*self.N:3*self.N, :]


        # Compute observables   
        self.population = self.C0 * np.sum(np.sqrt(np.arange(1, self.N+1))[:, np.newaxis] * self.f_n, axis=0)
        self.polarization = 0.5 * self.C0 * np.sum(np.sqrt(np.arange(1, self.N+1))[:, np.newaxis] * self.p_n, axis=0)
        
        # Compute absorption spectrum
        # self.T2_t = self.compute_T2(self.population)
        
        # Compute Fourier transform for absorption spectrum
        self.compute_absorption_spectrum()
    


    def compute_absorption_spectrum(self):
        """
        Compute the absorption spectrum from the time-dependent polarization.
        Absorption is proportional to Im[χ(ω)] = Im[P(ω)/E(ω)].
        """
        from scipy.fft import fft, fftfreq
        
        # Use FFT with proper normalization
        N = len(self.t)
        
        # Compute frequency grid
        freqs = fftfreq(N, self.dt)  # in 1/fs
        # Convert to energy in meV: E = ℏω = ℏ(2πf)
        energies = self.hbar * 2.0 * np.pi * freqs
        
        # Only keep positive frequencies
        pos_mask = (energies > 0) & (energies < 500)  # Focus on 0-150 meV
        
        # Apply window to reduce spectral leakage
        window = np.blackman(N)
        
        # Compute FFT with window
        P_fft = fft(self.polarization)
        E_array = np.array([self.laser_pulse(ti) for ti in self.t])
        E_fft = fft(E_array)
        
        # Extract positive frequencies
        P_fft_pos = P_fft[pos_mask]
        E_fft_pos = E_fft[pos_mask]
        energies_pos = energies[pos_mask]
        
        # Safe division with threshold
        E_max = np.max(np.abs(E_fft_pos))
        threshold = 1e-4 * E_max
        safe = np.abs(E_fft_pos) > threshold
        
        # Compute susceptibility and absorption
        chi = np.zeros_like(P_fft_pos, dtype=complex)
        chi[safe] = P_fft_pos[safe] / E_fft_pos[safe]
        
        # Absorption: Im[χ(ω)]
        alpha = np.imag(chi)
        
        self.spectrum_energy = energies_pos
        self.absorption_spectrum = alpha

        
 

    def extract_results(self):
        """
        Extract computed results after fitting.
        """
        return (self.t, self.f_n, self.p_n, self.population, self.polarization, self.absorption_spectrum, self.spectrum_energy)


    

