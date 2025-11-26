import numpy as np

def save_results(solver):
        np.savetxt(f'data/fn_chi{solver.chi_0}.dat', np.column_stack((solver.t[:, np.newaxis], solver.f_n)), delimiter='\t')
        np.savetxt(f'data/pn_chi{solver.chi_0}.dat', np.column_stack((solver.t[:, np.newaxis], solver.p_n.real, solver.p_n.imag)), delimiter='\t')
        np.savetxt(f'data/population_chi{solver.chi_0}.dat', np.column_stack((solver.t, solver.population)), fmt="%.2f", delimiter='\t')
        np.savetxt(f'data/polarization_chi{solver.chi_0}.dat', np.column_stack((solver.t, solver.polarization.real, solver.polarization.imag)), delimiter='\t')
        np.savetxt(f'data/absorption_spectrum_chi{solver.chi_0}.dat', np.column_stack((solver.spectrum_energy, solver.absorption_spectrum)), delimiter='\t')
def save_absorption_comparison(solver1, solver2):
    np.savetxt(f'data/absorption_comparison_with_chi{solver1.chi_0}.dat', np.column_stack((solver1.spectrum_energy, solver1.absorption_spectrum)), delimiter='\t')
    np.savetxt(f'data/absorption_comparison_without_chi{solver2.chi_0}.dat', np.column_stack((solver2.spectrum_energy, solver2.absorption_spectrum)), delimiter='\t')