import numpy as np
from numpy import pi, sqrt, log, exp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from tqdm import tqdm
import datetime
# Định nghĩa các hằng số
S = 1000
N = 1000                     # Số phương trình
hbar = 658.5                # meV.fs
E_R = 4.2                   # meV
delta_t = 10                # femtosecond
delta_0 = 20                # meV (Độ bơm để electron đi sâu vào dải -> delta_0 càng lớn thì electron càng vào sâu bên trong dải dẫn!) = hbar*omega_0 - E_g
t_max = 1000                 # femtosecond
t_0 = -3 * delta_t          # femtosecond
e_max = 300                 # meV
delta_e = e_max / N         # meV
dt = 2                      # femtosecond
khi_o = [0.001, 0.1, 0.2, 0.5]   # Tham số cường độ xung (0.1 - 2)
a0 = 125*1e-8               # Bán kính Bohr (cm)
gamma = 6.5*1e-20           # cm^{3}.fs^{-1}
T20 = 210
# Tính hệ số
C0 = (delta_e * sqrt(delta_e)) / (2 * (pi**2) * E_R**(3/2) * a0**3)

time_steps = np.arange(t_0, t_max, dt)
print(len(time_steps))
omega = np.linspace(-200, 200, S)
Et = np.zeros(len(time_steps))
Re_energy = np.zeros(S)
energy = np.zeros(S, dtype='complex')
p_omega = np.zeros(S, dtype='complex')
a_omega = np.zeros(S)

# Hàm tính toán
sqrt_n = np.sqrt(np.arange(1, N + 1))
g = np.zeros((N, N))    
for n in range(N):
    for n1 in range(N):
        if n1 != n:
            g[n, n1] = (1 / sqrt((n + 1) * delta_e)) * log(abs(((sqrt_n[n]) + sqrt_n[n1]) / ((sqrt_n[n] - sqrt_n[n1]))))

def En(Y, n): # Năng lượng rời rạc
    n1 = np.arange(N) != n  
    tong = np.dot(g[n, n1], ((Y[0][n1]).real + (Y[0][n1]).imag))
    return (sqrt(E_R) / pi) * delta_e * tong

def Omega(Y, n, t, khi0): # Tần số tái chuẩn hóa Rabi
    n1 = np.arange(N) != n  
    tong = np.dot(g[n, n1], Y[1][n1])
    Omega1 = (sqrt(pi) / (2 * delta_t)) * khi0 * exp(-(t / delta_t) ** 2)
    Omega2 = (sqrt(E_R) / (hbar * pi)) * delta_e * tong
    return (Omega1 + Omega2)

def Ft(t, Y, T2, khi0):
    F = np.zeros((2, N), dtype='complex')
    for i in range(N):
        F[0][i] = -2 * (Omega(Y, i, t, khi0) * np.conj(Y[1][i])).imag + 1j * -2 * (Omega(Y, i, t, khi0) * np.conj(Y[1][i])).imag
        F[1][i] = (-1j / hbar) * ((i + 1) * delta_e - delta_0 - En(Y, i)) * (Y[1][i]) \
                  + 1j * (1 - np.real(Y[0][i]) - np.imag(Y[0][i])) * Omega(Y, i, t, khi0) \
                  - (Y[1][i] / T2)
    return F

def RK4(Ft, t, Y, khi0, T2):
    Nt = C0 *  np.sum(sqrt_n[:N] * Y[0].real)
    T2 = 1 / (1 / T20 + gamma * Nt)

    k1 = dt * Ft(t, Y, T2, khi0)
    k2 = dt * Ft(t + dt / 2, Y + k1 / 2, T2, khi0)
    k3 = dt * Ft(t + dt / 2, Y + k2 / 2, T2, khi0)
    k4 = dt * Ft(t + dt, Y + k3, T2, khi0)

    Y = Y + (k1 + 2*k2 + 2*k3 +k4)/6

    Pt = delta_e * sqrt(delta_e) * np.sum(sqrt_n[:N] * (Y[1]))
    PT = delta_e * sqrt(delta_e) * np.sum(sqrt_n[:N] * np.abs(Y[1][0]))

    return Y, Nt, T2, Pt, PT

start = time.time()
##################################### FOURIER TRANSFORMATION #####################################
def Fourier(time_steps, Pt, khi0, file_name):
    
    for o in range(S):
        summ1 = 0 + 0j
        summ2 = 0 + 0j
        for t in range(len(time_steps)):
            summ1 += dt * khi0 * exp(-(time_steps[t]*time_steps[t]) / (delta_t*delta_t)) * np.exp(1j * omega[o] * time_steps[t] / hbar)
            summ2 += dt * Pt[t] * np.exp(1j * omega[o] * time_steps[t] / hbar)
           
        energy[o] = summ1
        Re_energy[o] = np.real(summ1)
        p_omega[o] = summ2
        a_omega[o] = np.imag(p_omega[o] / energy[o])

        write_FT(file_name, omega, energy, p_omega, a_omega)

    return a_omega

def write_result(filename, step, epsilon, Y1_real, Y2_abs):
    with open(filename, 'w') as file:
        file.write(f"{'#Time t':^20} {'Epsilon':^20} {'Y1_real':^20} {'Y2_abs':^40} \n")
        for i, t in enumerate(time_steps): 
            for j in range(len(Y1_real[i])):  
                file.write(f"{t:^20.5e} {epsilon[j]:^20.5e} {Y1_real[i][j]:^20.5e} {Y2_abs[i][j]:^40.5e} \n")
            file.write("\n")

def write_Nt_Pt(filename, step, Nt, Pt, PT, T2):
    with open(filename, 'w') as file:
        file.write(f"{'Time t':^20} {'Nt':^20} {'Pt (Real + Imag)':^30} {'|PT|':^20} {'T2':^20}\n")
        for t in range(len(time_steps)):
            file.write(f"{step[t]:^20.5e} "  
                       f"{Nt[t]:^20.5e} "       
                       f"{Pt[t].real:^15.5e} + {Pt[t].imag:^15.5e}j " 
                       f"{PT[t]:^20.5e} "        
                       f"{T2[t]:^20.5e}\n")   

def write_FT(filename, w, E, P, alpha):
    with open(filename, 'w') as file:
        file.write(f"#{'omega':^20} {'E':^20} {'P':^30} {'alpha':^30}\n ")
        for i in range(len(omega)):
            file.write(f"{w[i]:^20.5e} {E[i]:^20e} {P[i]:^30e} {alpha[i]:^30e} \n")


a_omega_collection = []
for khi in (khi_o):
    Y1_e_collection1 = []
    Y2_p_collection1 = []
    total_density_arr = []
    total_polarz1_arr = []
    abs_total_polarz1_arr = []
    T2_arr = []
    Nt = np.zeros(len(time_steps))
    Pt = np.zeros(len(time_steps), dtype = 'complex')
    PT = np.zeros(len(time_steps))
    T2_values = np.zeros(len(time_steps))

    epsilon = [(i + 1) * delta_e for i in range(N)] 

    Y = np.zeros((2, N), dtype='complex')
    Y[0][0] = 0 + 1j * 0  
    Y[1][0] = 0 + 1j * 0  
    
    for t in tqdm(range(len(time_steps))):
        Y, Nt[t], T2_values[t], Pt[t], PT[t] = RK4(Ft, time_steps[t], Y, khi, T20)

        Y1_e = []
        Y2_p = []

        for i in range(N):
            Y1_e.append(float(np.real(Y[0][i]))) 
            Y2_p.append(abs(Y[1][i]))  
        

        Y1_e_collection1.append(Y1_e)
        Y2_p_collection1.append(Y2_p)
        
    #write_result(f'SBE_khi0={khi}_T20=210.txt', time_steps, epsilon, Y1_e_collection1, Y2_p_collection1)
    #write_Nt_Pt(f'Nt Pt cung voi alpha khi0 = {khi} T20 = 210.txt', time_steps, Nt, np.abs(Pt), PT, T2_values)

    time_now = datetime.datetime.now().strftime("%H%M")
    file_name = f"alpha(omega)_khi0_{khi}_deltat_{delta_t}_delta0_{delta_0}_T2_{T20}_time_{time_now}.txt"
    a_omega = Fourier(time_steps, Pt, khi, file_name)
    a_omega_collection.append(a_omega)

plt.figure(figsize=(10, 6))
styles = ['-', '--', '=', '+']  
for i, khi in enumerate(khi_o):
    plt.plot(omega, a_omega_collection[i], label=f'khi0 = {khi}', linestyle=styles[i % len(styles)])

# Customize the plot
plt.title('Absorption Spectra: a_omega vs omega for different khi0')
plt.xlabel('omega (meV)')
plt.ylabel('a_omega (Imaginary part of p_omega / E)')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()


end = time.time()
interval = end - start
print("Quá trình chạy kết thúc trong ", interval)