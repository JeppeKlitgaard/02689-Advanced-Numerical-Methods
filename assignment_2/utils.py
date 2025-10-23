from math import gamma

import numpy as np
import numpy.typing as npt
from scipy import linalg
from scipy.special import eval_jacobi


def JacobiGQ(alpha: float, beta: float, N: int):
    """
    Compute N'th order Gauss-Jacobi quadrature nodes x and weights w
    for weight (1-x)^alpha (1+x)^beta on [-1, 1], with alpha,beta > -1.
    """
    if N == 0:
        x = np.array([-(alpha - beta) / (alpha + beta + 2.0)], dtype=float)
        w = np.array([2.0], dtype=float)  # matches the provided MATLAB code
        return x, w

    J = np.zeros((N + 1, N + 1), dtype=float)
    k = np.arange(0, N + 1, dtype=float)
    h1 = 2.0 * k + alpha + beta

    J[np.arange(N + 1), np.arange(N + 1)] = -0.5 * (alpha**2 - beta**2) / ((h1 + 2.0) * h1)
    if (alpha + beta) < 10.0 * np.finfo(float).eps:
        J[0, 0] = 0.0  # Legendre limit

    kk = np.arange(1, N + 1, dtype=float)
    b = 2.0 / (h1[:-1] + 2.0) * np.sqrt(
        kk * (kk + alpha + beta) * (kk + alpha) * (kk + beta)
        / ((h1[:-1] + 1.0) * (h1[:-1] + 3.0))
    )
    J += np.diag(b, 1)
    J += np.diag(b, -1)
    evals, evecs = np.linalg.eigh(J)
    x = evals  # Gauss nodes
    # Weights from first row of normalized eigenvectors
    # mu0 = integral_{-1}^{1} (1-x)^alpha (1+x)^beta dx
    mu0 = (2.0 ** (alpha + beta + 1.0)) * gamma(alpha + 1.0) * gamma(beta + 1.0) / gamma(alpha + beta + 2.0)
    w = mu0 * (evecs[0, :] ** 2)

    idx = np.argsort(x)
    x, w = x[idx], w[idx]

    return x, w


def JacobiGL(alpha: float, beta: float, N: int) -> npt.NDArray:
    """
    Returns the Gauss-Lobatto nodes, x ∈ [-1,1].

    You may need to do a change of basis to get a suitable interval
    for your problem.
    """
    x = np.zeros(N+1, dtype=float)
    if N == 1:
        x = np.array([-1, 1])
        return x
    xint, _ = JacobiGQ(alpha + 1, beta + 1, N-2)
    x = np.concatenate((np.array([-1]), xint, np.array([1])))
    return x


def JacobiP(x: np.ndarray, alpha: float, beta: float, N: int) -> np.ndarray:
    pass


def GradJacobiP(x: np.ndarray, alpha: float, beta: float, N: int) -> np.ndarray:
    pass


def constructV(x: np.ndarray, alpha: float = 0.0, beta: float = 0.0):
    "Function for constructing the Vandermonde matrix, V_{ij} = P_j(x_i)"
    N = x.shape[0]
    V = np.zeros((N, N))
    for n in range(N):
        V[n, :] = eval_jacobi(n, alpha, beta, x)
    return V.T


def constructVx(x: np.ndarray, alpha: float = 0.0, beta: float = 0.0):
    "Function for constructing the derivative of Vandermonde matrix, Vx = dV/dx"
    N = x.shape[0]
    Vx = np.zeros((N, N))
    for n in range(1, N):
        Vx[n, :] = ((n + alpha + beta + 1)/2) * eval_jacobi(n-1, alpha+1, beta+1, x)
    return Vx.T


def constructD(x: npt.NDArray, alpha: float = 0.0, beta: float = 0.0) -> npt.NDArray:
    "Function for constructing the differentiation matrix D = Vx V^{-1}"
    V = constructV(x, alpha=alpha, beta=beta)
    Vx = constructVx(x, alpha=alpha, beta=beta)
    # solve D = Vx V^{-1} <=> D V = Vx <=> V.T D.T = Vx.T <=> D.T = solve(V.T, Vx.T) <=> D = solve(V.T, Vx.T).T
    D = linalg.solve(V.T, Vx.T).T
    return D


def z2x(z, x0, xN): return ((xN - x0) * z + (xN + x0)) / 2
def x2z(x, x0, xN): return (2 * x - (xN + x0)) / (xN - x0)

def z2x_dx(x0, xN): return (xN - x0) / 2
def x2z_dz(x0, xN): return 2 / (xN - x0)