import numpy as np

def save_results(solver):
        np.savetxt('data/fn.dat', np.column_stack((solver.t, solver.f_n.T)), delimiter='\t')
        np.savetxt('data/pn.dat', np.column_stack((solver.t, solver.p_n.real.T, solver.p_n.imag.T)), delimiter='\t')
        np.savetxt('data/population.dat', np.column_stack((solver.t, solver.population)), fmt="%.2f", delimiter='\t')
        np.savetxt('data/polarization.dat', np.column_stack((solver.t, solver.polarization.real, solver.polarization.imag)), delimiter='\t')
        np.savetxt('data/absorption_spectrum.dat', np.column_stack((solver.spectrum_energy, solver.absorption_spectrum)), delimiter='\t')


def save_absorption_comparison(solver1, solver2):
    np.savetxt('data/absorption_comparison_with.dat', np.column_stack((solver1.spectrum_energy, solver1.absorption_spectrum)), delimiter='\t')
    np.savetxt('data/absorption_comparison_without.dat', np.column_stack((solver2.spectrum_energy, solver2.absorption_spectrum)), delimiter='\t')