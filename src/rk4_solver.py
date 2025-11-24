"""
RK4 (Runge-Kutta 4th order) ODE solver.
"""

import numpy as np


def rk4_solve(func, t_span, y0, t_eval):
    """
    Solve ODE using RK4 method.
    
    dy/dt = func(t, y)
    """
    n_points = len(t_eval)
    n_vars = len(y0)
    y_solution = np.zeros((n_vars, n_points))
    y_solution[:, 0] = y0
    
    t = t_eval[0]
    y = y0.copy()
    
    for i in range(1, n_points):
        dt = t_eval[i] - t_eval[i-1]
        
        k1 = func(t, y)
        k2 = func(t + 0.5*dt, y + 0.5*dt*k1)
        k3 = func(t + 0.5*dt, y + 0.5*dt*k2)
        k4 = func(t + dt, y + dt*k3)
        
        y = y + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        y_solution[:, i] = y
        t = t_eval[i]
    
    return {
        't': t_eval,
        'y': y_solution,
    }
