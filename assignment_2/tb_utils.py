import numpy as np
from math import gamma
from scipy.special import jacobi, legendre, roots_jacobi, eval_jacobi, factorial, gamma
from scipy import linalg
import numpy.typing as npt

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


def JacobiGL(alpha, beta, N):
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


def z2x(z, a, b): return ((b - a) * z + (b + a)) / 2
def x2z(x, a, b): return (2 * x - (b + a)) / (b - a)

def z2x_dz(a, b): return (b - a) / 2
def x2z_dx(a, b): return 2 / (b - a)


def constructV(z: np.ndarray, alpha=0, beta=0, a = -1, b = 1):
    "Function for constructing the Vandermonde matrix, V_{ij} = P_j(x_i)"
    N = z.shape[0]
    V = np.zeros((N, N))
    for n in range(N):
        V[n, :] = eval_jacobi(n, alpha, beta, z)
    # V = np.sqrt(2 / (b - a)) * V
    # V = V * np.sqrt((2*np.arange(1,N+1))/(b - a))
    return V.T


def constructVx(z: np.ndarray, alpha=0, beta=0, a = -1, b = 1):
    "Function for constructing the derivative of Vandermonde matrix, Vx = dV/dx; chain rule applied"
    N = z.shape[0]
    Vx = np.zeros((N, N))
    for n in range(1, N):
        Vx[n, :] = ((n + alpha + beta + 1)/2) * eval_jacobi(n-1, alpha+1, beta+1, z)
    Vx = 2 / (b - a) * Vx
    return Vx.T


def constructD(z: np.ndarray, alpha = 0, beta = 0, a = -1, b = 1, V = None, Vx = None):
    "Function for constructing the differentiation matrix D = Vx V^{-1}"
    if V is None:
        V = constructV(z=z, alpha=alpha, beta=beta, a=a, b=b)
    if Vx is None:
        Vx = constructVx(z=z, alpha=alpha, beta=beta, a=a, b=b)
    # solve D = Vx V^{-1} <=> D V = Vx <=> V.T D.T = Vx.T <=> D.T = solve(V.T, Vx.T) <=> D = solve(V.T, Vx.T).T
    D = linalg.solve(V.T, Vx.T).T
    return D


def jacobi_p_normalisation_const(alpha: float, beta: float, n: int):
    "Author: Jeppe"
    return (
        2 ** (alpha + beta + 1)
        * (gamma(n + alpha + 1) * gamma(n + beta + 1))
        / (factorial(n) * (2 * n + alpha + beta + 1) * gamma(n + alpha + beta + 1))
    )


def constructMinv(z = None, V = None, alpha=0, beta=0, a = -1, b = 1):
    "Function for constructing the mass matrix inverse M^{-1}"
    if V is None:
        V = constructV(z=z, alpha=alpha, beta=beta, a=a, b=b)
    N = V.shape[0]
    norm_const = jacobi_p_normalisation_const(alpha, beta, np.arange(N))
    V_normalized = V / np.sqrt(norm_const)
    M_inv = V_normalized @ V_normalized.T
    M_inv = M_inv * (2 / (b - a))
    return M_inv


def evaluate_jacobi_grid(x_eval, x_nodes, alpha=0, beta=0, a=-1, b=1):
    M = x_eval.shape[0]
    N = x_nodes.shape[0]
    A = np.zeros((N, M))
    z_eval = x2z(x_eval, a, b)
    for n in range(N):
        # A[n, :] = z2x(eval_jacobi(n, alpha, beta, z_eval), x0, xN)
        A[n, :] = eval_jacobi(n, alpha, beta, z_eval)
    return A.T
