import numpy as np
from math import gamma

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