#%%
import numpy as np
from numpy import pi, exp, sin, cos, sqrt, real, imag, conjugate
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib import rcParams
import matplotlib.transforms as mtransforms
import sys
import multiprocessing as mp

def on_press(event):
    print("my position:" ,event.button,event.xdata, event.ydata)

fontsize = 9
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['figure.dpi'] = 150

# serif:
plt.rcParams['pdf.fonttype'] = 42
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath, txfonts}') # txfonts
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'

# ================================================================
# version 3
# Solve the SPE equivalent to the Burgers equation using explicit/implicit FDM.
# The spatial derivative terms are approximated by second-order central difference.
# Time integral: first-order Euler method or third-order Runge-Kutta scheme.
# h = 1, m = 1
# update log:
# [231228][v3] a new set of \vec{f} and V_r
# [231230][v3] add third-order Runge-Kutta scheme to time integrate
# ================================================================

def dfdx(data, dx, nx):
    diff = np.zeros(nx, dtype='complex_')
    # Second order central difference:
    for i in range(1, nx-1):
        diff[i] = (data[i+1] - data[i-1]) / (2*dx)
    diff[0] = (data[1] - data[nx-1]) / (2*dx)
    diff[nx-1] = (data[0] - data[nx-2]) / (2*dx)
    return diff

def d2fdx2(data, dx, nx):
    diff = np.zeros(nx, dtype='complex_')
    # FDM:
    for i in range(1, nx-1):
        diff[i] = (data[i+1] - 2*data[i] + data[i-1]) / dx**2
    diff[0] = (data[nx-1] - 2*data[0] + data[1]) / dx**2
    diff[nx-1] = (data[nx-2] - 2*data[nx-1] + data[0]) / dx**2
    return diff

def compute_s(psi1, psi2):
    s1 = real(abs(psi1)**2 - abs(psi2)**2)
    s2 = real(1j * (conjugate(psi1)*psi2 - psi1*conjugate(psi2)))
    s3 = real(conjugate(psi1)*psi2 + psi1*conjugate(psi2))
    return s1, s2, s3

def compute_m(psi1, psi2, dx, nx):
    m1 = abs(dfdx(psi1, dx, nx))**2 - abs(dfdx(psi2, dx, nx))**2
    m2 = 1j * (conjugate(dfdx(psi1, dx, nx))*dfdx(psi2, dx, nx) - dfdx(psi1, dx, nx)*conjugate(dfdx(psi2, dx, nx)))
    m3 = conjugate(dfdx(psi1, dx, nx))*dfdx(psi2, dx, nx) - dfdx(psi1, dx, nx)*conjugate(dfdx(psi2, dx, nx))
    return m1, m2, m3

def compute_rho_velocity(psi1, psi2, dx, nx):
    rho = np.abs(psi1)**2 + np.abs(psi2)**2
    u = real((real(psi1)*dfdx(imag(psi1), dx, nx) - imag(psi1)*dfdx(real(psi1), dx, nx) + real(psi2)*dfdx(imag(psi2), dx, nx) - imag(psi2)*dfdx(real(psi2), dx, nx))) / rho
    return rho, u

def compute_f(psi1, psi2, s1, s2, s3, rho, u, dx, nx, nu):
    grads2 = dfdx(s1, dx, nx)**2 + dfdx(s2, dx, nx)**2 + dfdx(s3, dx, nx)**2
    tmp1 = -1/(4*rho)*(dfdx(rho, dx, nx)**2 - 2*rho*d2fdx2(rho, dx, nx) + grads2)
    tmp2 = 2*nu*(dfdx(rho, dx, nx)*dfdx(u, dx, nx) + (abs(dfdx(psi1, dx, nx))**2 + abs(dfdx(psi2, dx, nx))**2)*u)
    lam1 = 1/(rho**2*(grads2 - dfdx(rho, dx, nx)**2)) * (grads2*tmp1 - rho*dfdx(rho, dx, nx)*tmp2)
    lam2 = 1/(rho**2*(grads2 - dfdx(rho, dx, nx)**2)) * (-rho*dfdx(rho, dx, nx)*tmp1 + rho**2*tmp2)
    f1 = lam1 * s1 + lam2 * dfdx(s1, dx, nx)
    f2 = lam1 * s2 + lam2 * dfdx(s2, dx, nx)
    f3 = lam1 * s3 + lam2 * dfdx(s3, dx, nx)
    return f1, f2, f3

def compute_potential(psi1, psi2, dx, nx, nu):
    rho, u = compute_rho_velocity(psi1, psi2, dx, nx)
    s1, s2, s3 = compute_s(psi1, psi2)
    f1, f2, f3 = compute_f(psi1, psi2, s1, s2, s3, rho, u, dx, nx, nu)

    grads2 = dfdx(s1, dx, nx)**2 + dfdx(s2, dx, nx)**2 + dfdx(s3, dx, nx)**2
    Vr = -1/(4*rho**2)*(dfdx(rho, dx, nx)**2 - 2*rho*d2fdx2(rho, dx, nx) + grads2/2)
    Vi = nu*(2*(abs(dfdx(psi1, dx, nx))**2 + abs(dfdx(psi2, dx, nx))**2) - d2fdx2(rho, dx, nx)) / (2*rho)
    P = d2fdx2(s1, dx, nx)/(4*rho) - f1
    Q = (d2fdx2(s3, dx, nx) + 1j*d2fdx2(s2, dx, nx)) / (4*rho) - (f3 + 1j*f2)
    return Vr, Vi, P, Q

def evolution(psi1, psi2, dx, dt, nx, nu):
    # First Euler scheme:
    Vr, Vi, P, Q = compute_potential(psi1, psi2, dx, nx, nu)
    source = 1j*dt*((Vr + P + 1j*Vi)*psi1 + Q*psi2) - psi1
    psi1_new = LUD_LinearSovler(nx, a, b, c, source)
    source = 1j*dt*((Vr - P + 1j*Vi)*psi2 + conjugate(Q)*psi1) - psi2
    psi2_new = LUD_LinearSovler(nx, a, b, c, source)

    psi1_new = psi1_new / np.sqrt(np.sum(abs(psi1)**2 + abs(psi2)**2)) * np.sqrt(nx)
    psi2_new = psi2_new / np.sqrt(np.sum(abs(psi1)**2 + abs(psi2)**2)) * np.sqrt(nx)
    return psi1_new, psi2_new

from scipy.signal import savgol_filter
def smoothing(data):
    data_filter = savgol_filter(real(data), 11, 5) + 1j*savgol_filter(imag(data), 11, 5)
    return data_filter

def LUD_LinearSovler(N, a, b, c, y):
    p = np.zeros(N, dtype='complex_')
    q = np.zeros(N-2, dtype='complex_')
    r = np.zeros(N-1, dtype='complex_')
    p[0] = b[0]
    q[0] = c[0] / p[0]
    r[0] = a[0] / p[0]
    for i in range(1, N-2):
        p[i] = b[i] - a[i]*q[i-1]
        q[i] = c[i] / p[i]
        r[i] = -a[i]*r[i-1] / p[i]
    p[N-2] = b[N-2] - a[N-2]*q[N-3]
    r[N-2] = (c[N-2] - a[N-2]*r[N-3]) / p[N-2]
    p[N-1] = b[N-1] - a[N-1]*r[N-2]
    s = c[N-1] / p[N-1]

    t = np.zeros(N-1, dtype='complex_')
    t[N-2] = r[N-2]
    for i in range(N-3, -1, -1):
        t[i] = r[i] - q[i]*t[i+1]

    v = np.zeros(N, dtype='complex_')
    v[0] = y[0] / p[0]
    for i in range(1, N):
        v[i] = (y[i] - a[i]*v[i-1]) / p[i]

    w = np.zeros(N, dtype='complex_')
    w[N-1] = v[N-1]
    w[N-2] = v[N-2]
    for i in range(N-3, -1, -1):
        w[i] = v[i] - q[i]*w[i+1]

    sol = np.zeros(N, dtype='complex_')
    sol[N-1] = (w[N-1] - w[0]*s) / (1 - t[0]*s)
    for i in range(N-1):
        sol[i] = w[i] - t[i]*sol[N-1]

    return sol

def test(psi1, psi2, dx, nx, nu):
    rho, u = compute_rho_velocity(psi1, psi2, dx, nx)
    s1, s2, s3 = compute_s(psi1, psi2)
    f1, f2, f3 = compute_f(psi1, psi2, s1, s2, s3, rho, u, dx, nx, nu)

    grads2 = dfdx(s1, dx, nx)**2 + dfdx(s2, dx, nx)**2 + dfdx(s3, dx, nx)**2
    Vr, Vi, P, Q = compute_potential(psi1, psi2, dx, nx, nu)
    tmp = 1/(2*rho)*(abs(dfdx(psi1, dx, nx))**2 + abs(dfdx(psi2, dx, nx))**2) - u**2/2 - Vr + (s1*f1 + s2*f2 + s3*f3)/rho
    rhs = dfdx(rho, dx, nx)/rho**2 * (s1*f1 + s2*f2 + s3*f3 + 1/(4*rho)*(dfdx(rho, dx, nx)**2 - 2*rho*d2fdx2(rho, dx, nx) + grads2 + 8*nu*rho**2*dfdx(u, dx, nx))) - (dfdx(s1, dx, nx)*f1 + dfdx(s2, dx, nx)*f2 + dfdx(s3, dx, nx)*f3)/rho + 2*nu/rho * (abs(dfdx(psi1, dx, nx))**2 + abs(dfdx(psi2, dx, nx))**2)*u + dfdx(tmp, dx, nx)
    print('The remaining terms of the equation:', np.max(abs(rhs)), np.max(abs(nu*d2fdx2(u, dx, nx))))
    return np.max(abs(rhs))

def monitor(psi1, psi2, dx, nx):
    rho, u = compute_rho_velocity(psi1, psi2, dx, nx)
    s1, s2, s3 = compute_s(psi1, psi2)
    f1, f2, f3 = compute_f(psi1, psi2, s1, s2, s3, rho, u, dx, nx, nu)

    if np.max(abs(rho)) >= 1e5:
        print('rho blow up')
        return 1
    elif np.max(abs(u)) >= 1e5:
        print('u blow up')
        return 1
    elif np.max(abs(s1)) >= 1e5 or np.max(abs(s2)) >= 1e5 or np.max(abs(s3)) >= 1e5:
        print('s blow up')
        return 1
    return 0


# set parameters:
nu = 1e-2
nx = 2**9
dx = 2*pi / (nx)
T = 1

x = np.arange(0, 2*pi, dx)

# initialize the wave function:
psi1 = np.zeros(nx, dtype='complex_')
psi2 = np.zeros(nx, dtype='complex_')
psi1 = cos(x) * (cos(cos(x)) - 1j*sin(cos(x)))
psi2 = sin(x) * (cos(cos(x)) - 1j*sin(cos(x)))

t = 0.0
dt = 0.001
dt_output = 0.01
t_output = 0
psi1_t_array = np.empty((0, nx+1))
psi2_t_array = np.empty((0, nx+1))
t_array = np.array([])
rhs_array = np.array([])
alpha = (nu + 1j/2)*dt/dx**2
a = alpha*np.ones(nx)
b = -(2*alpha + 1)*np.ones(nx)
c = a
while t < T + dt:
    if t >= t_output:
        print('===============================================')
        print('t =', t, '     dt =', dt, '     dt/dx =', dt/dx)
        t_output = t_output + dt_output
        tmp = psi1
        tmp = np.append(tmp, psi1[0])
        psi1_t_array = np.append(psi1_t_array, [tmp], axis=0)
        tmp = psi2
        tmp = np.append(tmp, psi2[0])
        psi2_t_array = np.append(psi2_t_array, [tmp], axis=0)
        t_array = np.append(t_array, t)
        rhs_max = test(psi1, psi2, dx, nx, nu)
        rhs_array = np.append(rhs_array, rhs_max)
        rho, u = compute_rho_velocity(psi1, psi2, dx, nx)
        print('total mass:', np.sum(rho))
    psi1, psi2 = evolution(psi1, psi2, dx, dt, nx, nu)
    psi1 = smoothing(psi1) # suppress numerical oscillations
    psi2 = smoothing(psi2)
    if monitor(psi1, psi2, dx, nx):
        break
    t = t + dt

U = np.zeros((np.size(psi1_t_array, 0), nx+1))
rho_t_array = np.zeros((np.size(psi1_t_array, 0), nx+1))
for i in range(np.size(psi1_t_array, 0)):
    rho_t_array[i, 0:-1], U[i, 0:-1] = compute_rho_velocity(psi1_t_array[i, 0:-1], psi2_t_array[i, 0:-1], dx, nx)
U[:, -1] = U[:, 0]
rho_t_array[:, -1] = rho_t_array[:, 0]
x = np.append(x, 2*pi)