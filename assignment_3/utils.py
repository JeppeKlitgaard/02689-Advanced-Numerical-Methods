import numpy as np
import numpy.typing as npt
from scipy.special import gamma


# Problem Functions
def get_hill_center(
    t: float, *, x_0: npt.NDArray, y_0: npt.NDArray
) -> tuple[npt.NDArray, npt.NDArray]:
    omega = 2 * np.pi
    x_t = x_0 * np.cos(omega * t) + y_0 * np.sin(omega * t)
    y_t = -x_0 * np.sin(omega * t) + y_0 * np.cos(omega * t)
    return x_t, y_t


def u_exact_hill(x, y, t, *, x_0, y_0, sigma):
    x_center_t, y_center_t = get_hill_center(t, x_0=x_0, y_0=y_0)
    return np.exp(-((x - x_center_t) ** 2 + (y - y_center_t) ** 2) / (2 * sigma**2))


def u_exact_hill_dt(x, y, t, *, x_0, y_0, sigma):
    omega = 2 * np.pi
    xc, yc = get_hill_center(t, x_0=x_0, y_0=y_0)
    xc_t = omega * yc
    yc_t = -omega * xc
    u = u_exact_hill(x, y, t, sigma=sigma, x_0=x_0, y_0=y_0)
    return u * ((x - xc) * xc_t + (y - yc) * yc_t) / (sigma**2)


def f_rhs_hill(x: npt.NDArray, y: npt.NDArray, t: float) -> npt.NDArray:
    return np.zeros_like(x)


def advection_velocity_field(x: npt.NDArray, y: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray]:
    return (2 * np.pi * y, -2 * np.pi * x)


def get_boundary_normal(x, y):
    norm = np.sqrt(x**2 + y**2)
    return x / norm, y / norm


def jacobi_gauss_quadrature(
    N: int, alpha: float, beta: float
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Compute N'th order Gauss-Jacobi quadrature nodes x and weights w
    for weight (1-x)^alpha (1+x)^beta on [-1, 1], with alpha,beta > -1.

    Adapted from APEK MATLAB code.
    """
    if N == 0:
        x = np.array([-(alpha - beta) / (alpha + beta + 2.0)], dtype=float)
        w = np.array([2.0], dtype=float)  # matches the provided MATLAB code
        return x, w

    J = np.zeros((N + 1, N + 1), dtype=float)
    k = np.arange(0, N + 1, dtype=float)
    h1 = 2.0 * k + alpha + beta

    J[np.arange(N + 1), np.arange(N + 1)] = (
        -0.5 * (alpha**2 - beta**2) / ((h1 + 2.0) * h1)
    )
    if (alpha + beta) < 10.0 * np.finfo(float).eps:
        J[0, 0] = 0.0  # Legendre limit

    kk = np.arange(1, N + 1, dtype=float)
    b = (
        2.0
        / (h1[:-1] + 2.0)
        * np.sqrt(
            kk
            * (kk + alpha + beta)
            * (kk + alpha)
            * (kk + beta)
            / ((h1[:-1] + 1.0) * (h1[:-1] + 3.0))
        )
    )
    J += np.diag(b, 1)
    J += np.diag(b, -1)
    evals, evecs = np.linalg.eigh(J)
    x = evals  # Gauss nodes
    # Weights from first row of normalized eigenvectors
    # mu0 = integral_{-1}^{1} (1-x)^alpha (1+x)^beta dx
    mu0 = (
        (2.0 ** (alpha + beta + 1.0))
        * gamma(alpha + 1.0)
        * gamma(beta + 1.0)
        / gamma(alpha + beta + 2.0)
    )
    w = mu0 * (evecs[0, :] ** 2)

    idx = np.argsort(x)
    x, w = x[idx], w[idx]

    return x, w


def jacobi_gauss_lobatto(
    N: int, alpha: float, beta: float
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Returns the N+1 Gauss-Lobatto nodes for the specified Jacobi polynomial.

    Gauss-Lobatto nodes include the endpoints -1 and 1.
    """
    x = np.zeros(N + 1, dtype=float)
    if N == 1:
        x = np.array([-1, 1])
        w = np.array([1, 1])
        return x, w
    
    xint, wint = jacobi_gauss_quadrature(N - 2, alpha + 1, beta + 1)
    x = np.concatenate([[-1.0], xint, [1.0]])
    w_boundary = 2.0 / (N * (N + 1))
    w = np.concatenate([[w_boundary], wint, [w_boundary]])
    return x, w