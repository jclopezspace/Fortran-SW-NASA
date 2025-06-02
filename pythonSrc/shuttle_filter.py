# Python translation of Fortran Shuttle Project
# Author: Steve Winward (original Fortran)
# Converted to Python by Copilot

import sys
import math

# ASCII NASA Shuttle
SHUTTLE_ART = r"""
           _
          /^\
          |-|
          | |
          |N|
          |A|
          |S|
          |A|
         /| |\
        /_|_|_\
           W
"""

def title():
    print("\n==============================")
    print("   NASA Shuttle Filter App   ")
    print("==============================\n")
    print(SHUTTLE_ART)

def condition_num(gamma):
    return (1 + abs(gamma)) / (1 - abs(gamma))

def get_rhs(gamma, n, s_unfiltered):
    a = (1.0 + gamma) / 2.0
    a_over_2 = a / 2.0
    rhs = [0.0] * (n-2)
    rhs[0] = a * s_unfiltered[1] + a_over_2 * s_unfiltered[2] + (a_over_2 - gamma/2.0) * s_unfiltered[0]
    for i in range(1, n-3):
        rhs[i] = a_over_2 * s_unfiltered[i] + a * s_unfiltered[i+1] + a_over_2 * s_unfiltered[i+2]
    rhs[n-3] = a_over_2 * s_unfiltered[n-3] + a * s_unfiltered[n-2] + (a_over_2 - gamma/2.0) * s_unfiltered[n-1]
    return rhs

def initialize_tridiag_matrix(gamma, n):
    u = [0.0] * n
    d = [1.0] * n
    l = [0.0] * n
    u[0] = gamma / 2.0
    for i in range(1, n-1):
        l[i] = gamma / 2.0
        u[i] = gamma / 2.0
    l[n-1] = gamma / 2.0
    u[n-1] = 0.0
    return l, d, u

def norm_infinity(vec):
    return max(abs(x) for x in vec)

def tri_diag_matrix_times_vector(l, d, u, x):
    n = len(x)
    result = [0.0] * n
    result[0] = d[0]*x[0] + u[0]*x[1]
    for i in range(1, n-1):
        result[i] = l[i]*x[i-1] + d[i]*x[i] + u[i]*x[i+1]
    result[n-1] = l[n-1]*x[n-2] + d[n-1]*x[n-1]
    return result

def jacobi_iteration(l, d, u, b, x_init, tol, max_iter):
    n = len(x_init)
    x_new = [0.0] * n
    residual = [0.0] * n
    b_norm = norm_infinity(b)
    tol_prod = tol * b_norm
    temp_vec = tri_diag_matrix_times_vector(l, d, u, x_init)
    residual = [b[i] - temp_vec[i] for i in range(n)]
    for k in range(max_iter):
        res_norm = norm_infinity(residual)
        if res_norm < tol_prod:
            return k+1
        x_new[0] = b[0] - u[0]*x_init[1]
        for i in range(1, n-1):
            x_new[i] = b[i] - l[i]*x_init[i-1] - u[i]*x_init[i+1]
        x_new[n-1] = b[n-1] - l[n-1]*x_init[n-2]
        for i in range(n):
            residual[i] = x_new[i] - x_init[i]
        x_init = x_new[:]
    return max_iter

def succ_over_relaxation(l, d, u, b, x_init, w, n, max_iter=1000, tol=1e-8):
    x = x_init[:]
    for k in range(max_iter):
        x_old = x[:]
        x[0] = (1-w)*x[0] + w/d[0]*(b[0] - u[0]*x[1])
        for i in range(1, n-1):
            x[i] = (1-w)*x[i] + w/d[i]*(b[i] - l[i]*x[i-1] - u[i]*x[i+1])
        x[n-1] = (1-w)*x[n-1] + w/d[n-1]*(b[n-1] - l[n-1]*x[n-2])
        if norm_infinity([x[i]-x_old[i] for i in range(n)]) < tol:
            return k+1
    return max_iter

def read_in_file(filename):
    with open(filename, 'r') as f:
        data = [float(line.strip()) for line in f if line.strip()]
    return len(data), data

def main():
    title()
    filename = "FILTER_IN.DAT"
    n, s = read_in_file(filename)
    print(f"The file {filename} was successfully loaded with {n} data values\n")
    l, d, u = initialize_tridiag_matrix(0.0, n-2)
    b = [0.0]*(n-2)
    x_init = [0.0]*(n-2)
    gamma = float(input("Please input a value for gamma (-1, 1): "))
    print(f"You entered gamma = {gamma}\n")
    b = get_rhs(gamma, n, s)
    l, d, u = initialize_tridiag_matrix(gamma, n-2)
    print(f"{'(W, Iterations)':>15}")
    w_vals = []
    iters = []
    for i in range(1, 100):
        w = 1.0 + i*0.01
        num_iter = succ_over_relaxation(l, d, u, b, x_init, w, n-2)
        w_vals.append(w)
        iters.append(num_iter)
        print(f"{w:4.2f} {num_iter:6d}")
    print("\nASCII Scatter Plot (Iterations vs W):\n")
    ascii_scatter_plot(w_vals, iters)

def ascii_scatter_plot(x_vals, y_vals, width=60, height=20):
    min_x, max_x = min(x_vals), max(x_vals)
    min_y, max_y = min(y_vals), max(y_vals)
    plot = [[' ' for _ in range(width)] for _ in range(height)]
    for x, y in zip(x_vals, y_vals):
        px = int((x - min_x) / (max_x - min_x) * (width-1))
        py = int((y - min_y) / (max_y - min_y) * (height-1))
        py = height - 1 - py
        plot[py][px] = '*'
    for row in plot:
        print(''.join(row))
    print(f"{' ' * (width//2 - 5)}W axis")
    print(f"Iterations (Y axis)")

if __name__ == "__main__":
    main()
