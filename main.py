"""
Semiconductor Bloch Equations (SBE) Solver

This script solves the multi-level SBEs using equations (0.5a) and (0.5b).
"""


from src.sbe_solver import SBESolver
from src.params import PARAMS
from src.utils import save_results, save_absorption_comparison
from src.visualization import plot_requested_plots, plot_absorption_comparison


def main():
    """Main function to set up and solve the multi-level SBE for GaAs."""
    
    # Get parameters from params.py
    chi_0 = PARAMS['chi_0']
    delta_t = PARAMS['delta_t']
    Delta_0 = PARAMS['Delta_0']
    E_R = PARAMS['E_R']
    T2 = PARAMS['T2']                    # Constant dephasing time
    T2_0 = PARAMS['T2_0']                # Initial dephasing time
    gamma = PARAMS['gamma']              # Scattering coefficient
    hbar = PARAMS['hbar']
    N = PARAMS['N']
    Delta_epsilon = PARAMS['Delta_epsilon']
    t_0 = PARAMS['t_0']
    t_max = PARAMS['t_max']
    dt = PARAMS['dt']
    C_0 = PARAMS['C_0']
    
    # Create and fit solver
    solver1 = SBESolver()
    solver2 = SBESolver()
    
    
    solver1.fit(
        chi_0=chi_0,
        delta_t=delta_t,
        Delta_0=Delta_0,
        E_R=E_R,
        T2=T2,
        T2_0=T2_0,      
        gamma=gamma,
        hbar=hbar,
        N=N,
        Delta_epsilon=Delta_epsilon,
        t_0=t_0,
        t_max=t_max,
        dt=dt,
        C_0=C_0,
        with_coulomb=True
    )

    solver2.fit(
        chi_0=chi_0,
        delta_t=delta_t,
        Delta_0=Delta_0,
        E_R=E_R,
        T2=T2,
        T2_0=T2_0,      
        gamma=gamma,
        hbar=hbar,
        N=N,
        Delta_epsilon=Delta_epsilon,
        t_0=t_0,
        t_max=t_max,
        dt=dt,
        C_0=C_0,
        with_coulomb=False
    )
    

    # Save solutions
    save_results(solver1)
    save_absorption_comparison(solver1, solver2)


    # Plot results
    plot_requested_plots(solver1)
    plot_absorption_comparison(solver1, solver2)
    


if __name__ == "__main__":
    main()
