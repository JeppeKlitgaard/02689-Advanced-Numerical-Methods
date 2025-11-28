import numpy as np
from scipy.special import eval_jacobi, factorial, gamma

def jacobi_p_normalisation_const(alpha: float, beta: float, n: int):
    return (
        2 ** (alpha + beta + 1)
        * (gamma(n + alpha + 1) * gamma(n + beta + 1))
        / (factorial(n) * (2 * n + alpha + beta + 1) * gamma(n + alpha + beta + 1))
    )

def JacobiP_orthonormal(x, alpha, beta, N):
    norm_constant = np.sqrt(jacobi_p_normalisation_const(alpha, beta, N))
    return eval_jacobi(N, alpha, beta, x) / norm_constant


def GradJacobiP_orthonormal(x, alpha, beta, n):
    if n == 0:
        return np.zeros_like(x)
    return np.sqrt(n*(n + alpha + beta + 1)) * JacobiP_orthonormal(x, alpha+1, beta+1, n-1)


def rs2ab(r, s, tol=1e-12):
    a = np.empty_like(r, dtype=float)
    mask = np.abs(1.0 - s) > tol
    a[mask] = 2.0 * (1.0 + r[mask]) / (1.0 - s[mask]) - 1.0
    a[~mask] = -1.0
    b = s
    return a, b


def Simplex2DP_orthonormal(a, b, i, j):
    # orthonormal basis, not just orthogonal, as in the book
    h1 = JacobiP_orthonormal(a, 0, 0, i) # / np.sqrt(jacobi_p_normalisation_const(0, 0, i))
    h2 = JacobiP_orthonormal(b, 2*i + 1, 0, j) #  / np.sqrt(jacobi_p_normalisation_const(2*i + 1, 0, j))
    P = np.sqrt(2.0) * h1 * h2 * (1 - b)**i
    return P


def Vandermonde2D_orthonormal(N, r, s):
    a, b = rs2ab(r, s)
    num_nodes = r.shape[0]
    num_modes = (N + 1)*(N + 2)//2
    V = np.zeros((num_nodes, num_modes))
    sk = 0
    for i in range(N + 1):
        for j in range(N - i + 1):
            V[:, sk] = Simplex2DP_orthonormal(a, b, i, j)
            sk += 1
    # V = V[LocalReorder, :][:, LocalReorder]
    return V


def GradSimplex2DP(a, b, id, jd):

    fa  = JacobiP_orthonormal(a, 0, 0, id)
    gb  = JacobiP_orthonormal(b, 2*id + 1, 0, jd)
    dfa = GradJacobiP_orthonormal(a, 0, 0, id)
    dgb = GradJacobiP_orthonormal(b, 2*id + 1, 0, jd)

    # r-derivative part
    if id == 0:
        dmodedr = dfa * gb  # dfa is zero
    else:
        dmodedr = dfa * gb * (0.5*(1 - b))**(id - 1)

    # s-derivative first part
    dmodeds = dfa * (gb * (0.5*(1 + a)))
    if id > 0:
        dmodeds *= (0.5*(1 - b))**(id - 1)

    # tmp term
    tmp = dgb * (0.5*(1 - b))**id
    if id > 0:
        tmp = tmp - 0.5*id*gb*(0.5*(1 - b))**(id - 1)

    dmodeds = dmodeds + fa * tmp

    # normalization
    factor = 2**(id + 0.5)
    dmodedr *= factor
    dmodeds *= factor

    return dmodedr, dmodeds


def GradVandermonde2D_orthonormal(N, r, s):
    a, b = rs2ab(r, s)
    num_nodes = r.shape[0]
    num_modes = (N + 1)*(N + 2)//2
    V2Dr = np.zeros((num_nodes, num_modes))
    V2Ds = np.zeros((num_nodes, num_modes))
    sk = 0
    for i in range(N + 1):
        for j in range(N - i + 1):
            dr, ds = GradSimplex2DP(a, b, i, j)
            V2Dr[:, sk] = dr
            V2Ds[:, sk] = ds
            sk += 1
    # V2Dr = V2Dr[LocalReorder, :][:, LocalReorder]
    # V2Ds = V2Ds[LocalReorder, :][:, LocalReorder]
    return V2Dr, V2Ds


def Dmatrices2D(V, Vr, Vs):
    Dr = np.linalg.solve(V.T, Vr.T).T
    Ds = np.linalg.solve(V.T, Vs.T).T
    return Dr, Ds


def GeometricFactors2D(x, y, Dr, Ds):
    """
    Compute metric terms rx,sx,ry,sy and Jacobian J at nodal points.
    x,y: (Np,) physical coordinates
    Dr, Ds: (Np,Np) differentiation matrices wrt reference (r,s)
    Returns rx, sx, ry, sy, J (each (Np,))
    """
    xr = Dr @ x
    xs = Ds @ x
    yr = Dr @ y
    ys = Ds @ y
    J = xr * ys - xs * yr
    rx = ys / J
    sx = -yr / J
    ry = -xs / J
    sy = xr / J
    return rx, sx, ry, sy, J


def Dmatrices2D_xy(Dr, Ds, rx, sx, ry, sy):
    # same as diag(rx) @ Dr
    D_x = rx[:, None] * Dr + sx[:, None] * Ds
    D_y = ry[:, None] * Dr + sy[:, None] * Ds
    return D_x, D_y
