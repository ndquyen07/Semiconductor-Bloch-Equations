import matplotlib.pyplot as plt
import numpy as np



def plot_requested_plots(solver):
    """
    Plot the specific requested plots: f_n vs t, p_n vs t, population vs t,
    polarization vs t, |polarization| vs t, and α vs energy.
    """
    t = solver.t
    f_n = solver.f_n
    p_n = solver.p_n
    N = solver.N
    
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    
    # 1. f_n follow t (occupation for first 5 levels)
    ax1 = axes[0, 0]
    colors = plt.cm.tab10(np.linspace(0, 1, 5))
    for i in range(min(5, N)):
        ax1.plot(t, f_n[:, i], label=f'f_{i+1}', color=colors[i], linewidth=2)
    ax1.set_xlabel('Time (fs)')
    ax1.set_ylabel('f_n')
    ax1.set_title('f_n vs t (first 5 levels)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. p_n follow t (microscopic polarization for first few levels)
    ax2 = axes[0, 1]
    colors = plt.cm.viridis(np.linspace(0, 1, min(1, N)))
    for i in range(min(1, N)):
        ax2.plot(t, p_n[:, i].real, label=f'Re(p_{i+1})', color=colors[i], linestyle='-')
        ax2.plot(t, p_n[:, i].imag, label=f'Im(p_{i+1})', color=colors[i], linestyle='--')
    ax2.set_xlabel('Time (fs)')
    ax2.set_ylabel('p_n')
    ax2.set_title('p_n vs t')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. population follow t
    ax3 = axes[1, 0]
    ax3.plot(t, solver.population, 'b-', linewidth=2)
    ax3.set_xlabel('Time (fs)')
    ax3.set_ylabel('Population')
    ax3.set_title('Population vs t')
    ax3.grid(True, alpha=0.3)
    
    # 4. polarization follow t
    ax4 = axes[1, 1]
    ax4.plot(t, solver.polarization.real, 'b-', label='Re(P)', linewidth=2)
    ax4.plot(t, solver.polarization.imag, 'r--', label='Im(P)', linewidth=2)
    ax4.set_xlabel('Time (fs)')
    ax4.set_ylabel('Polarization')
    ax4.set_title('Polarization vs t')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. |polarization| follow t
    ax5 = axes[2, 0]
    pol_mag = np.abs(solver.polarization)
    ax5.plot(t, pol_mag, 'g-', linewidth=2)
    ax5.set_xlabel('Time (fs)')
    ax5.set_ylabel('|Polarization|')
    ax5.set_title('|Polarization| vs t')
    ax5.grid(True, alpha=0.3)
    
    # 6. α follow energy
    ax6 = axes[2, 1]
    ax6.plot(solver.spectrum_energy, solver.absorption_spectrum, 'b-', linewidth=2)
    ax6.set_xlabel('Energy (meV)')
    ax6.set_ylabel('α(ω)')
    ax6.set_title('α vs Energy')
    ax6.grid(True, alpha=0.3)
    ax6.legend()
    
    plt.tight_layout()
    plt.savefig('result/sbe.png', dpi=300, bbox_inches='tight')
    plt.show()


def plot_absorption_comparison(solver_with, solver_without):
    """
    Plot comparison of absorption spectra with and without Coulomb interactions.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(solver_with.spectrum_energy, solver_with.absorption_spectrum, 'b-', linewidth=2, label='With Coulomb')
    ax.plot(solver_without.spectrum_energy, solver_without.absorption_spectrum, 'r--', linewidth=2, label='Without Coulomb')
    
    ax.set_xlabel('Energy (meV)', fontsize=12)
    ax.set_ylabel('Absorption α(ω)', fontsize=12)
    ax.set_title('Absorption Spectrum: With vs Without Coulomb Interactions', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('result/absorption_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()