"""
Core solver for Semiconductor Bloch Equations.
"""

import numpy as np
from src.rk4_solver import rk4_solve

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

        E_n_array = (np.sqrt(self.E_R) / np.pi) * self.Delta_epsilon * (self.g_matrix @ (2.0 * f_levels))
        
        return E_n_array
    
        
    

    def laser_pulse(self, t):
        """Gaussian laser pulse.
        """
        return 0.5 * (np.sqrt(np.pi) / self.delta_t) * self.chi_0 * np.exp(-t**2 / self.delta_t**2)
    

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
        pulse = self.laser_pulse(t)
        coulomb_term = (1.0 / self.hbar) * (np.sqrt(self.E_R) / np.pi) * self.Delta_epsilon * self.g_matrix @ p_levels

        if self.with_coulomb:
            Omega_R_array = pulse + coulomb_term
        else:
            Omega_R_array = pulse
        
        return Omega_R_array
    

    def derivatives(self, t, y):
        """
        Compute the derivatives for the SBE system with N energy levels.
        """


        # Extract state variables
        f_n = y[:self.N]
        p_n = y[self.N : 2*self.N] + 1j * y[2*self.N : 3*self.N]


        # Compute current population and time-dependent T2
        population = self.C0 * np.sum(np.sqrt(np.arange(1, self.N+1)) * f_n)
        T2_current = self.compute_T2(population)

        if self.with_coulomb:
            E_n_array = self.coulomb_energy(f_n)
        else:
            E_n_array = np.zeros_like(f_n)

            
        Omega_R_array = self.rabi_frequency(t, p_n)

        # Compute derivatives
        dfdt = -2.0 * np.imag(Omega_R_array * np.conj(p_n))
        
        n_array = np.arange(1, self.N + 1)
        detuning = n_array * self.Delta_epsilon - self.Delta_0 - E_n_array
        dpdt = (-1j / self.hbar) * detuning * p_n + 1j * (1.0 - 2.0 * f_n) * Omega_R_array - p_n / T2_current
        
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

        # Initial conditions: all populations and polarizations zero
        y0 = np.zeros(3 * self.N)

        # Time grid
        t_points = int((self.t_max - self.t_0) / self.dt) + 1
        t_eval = np.linspace(self.t_0, self.t_max, t_points)
        
        # Solve the SBEs using RK4
        solution = rk4_solve(self.derivatives, y0, t_eval)
        
        # Extract results
        self.t = solution['t']
        self.f_n = solution['y'][:, :self.N]
        self.p_n = solution['y'][:, self.N:2*self.N] + 1j * solution['y'][:, 2*self.N:3*self.N]


        # Compute observables   
        self.population = self.C0 * np.sum(np.sqrt(np.arange(1, self.N+1)) * self.f_n, axis=1)
        self.polarization = 0.5 * self.C0 * np.sum(np.sqrt(np.arange(1, self.N+1)) * self.p_n, axis=1)

        # Compute Fourier transform for absorption spectrum
        self.compute_absorption_spectrum()
    


    def fourier_transform(self, signal):
        N = len(self.t)
        omega = np.linspace(-200, 200, N)
        
        omega_grid, t_grid = np.meshgrid(omega, self.t, indexing='ij')
        
        signal_fft = self.dt * np.sum(signal * np.exp(1j/self.hbar * omega_grid * t_grid), axis=1)
        
        return omega, signal_fft



    def compute_absorption_spectrum(self):
        omega, P_omega =  self.fourier_transform(self.polarization)
        _, E_omega = self.fourier_transform(self.laser_pulse(self.t))


        self.spectrum_energy = omega
        self.absorption_spectrum = np.imag(P_omega / E_omega)

    

    def extract_results(self):
        """
        Extract computed results after fitting.
        """
        return (self.t, self.f_n, self.p_n, self.population, self.polarization, self.absorption_spectrum, self.spectrum_energy)


    
if __name__ == "__main__":
    import numpy as np

    n = np.array([1, 3, 2])
    n1 = np.array([1, 2, 3])
    print(n.reshape(-1, 1) * n1)
    print(np.sum(n.reshape(-1, 1) * n1, axis=1))


    # print("n_grid:\n", n_grid)
    # print("n1_grid:\n", n1_grid)





