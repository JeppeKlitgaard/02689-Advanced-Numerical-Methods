import numpy as np
import numpy.typing as npt
from scipy.special import gamma, factorial
from typing import Protocol


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


def jacobi_gauss_lobatto(N: int, alpha: float, beta: float) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Returns the N+1 Gauss-Lobatto nodes for the specified Jacobi polynomial.

    Gauss-Lobatto nodes include the endpoints -1 and 1.
    """
    x = np.zeros(N + 1, dtype=float)
    if N == 1:
        x = np.array([-1, 1])
        w = np.array([1, 1])
        return x, w
    xint, wint = jacobi_gauss_quadrature(N-2, alpha + 1, beta + 1)
    x = np.concatenate([[-1.0], xint, [1.0]])
    w_boundary = 2.0 / (N * (N + 1))
    w = np.concatenate([[w_boundary], wint, [w_boundary]])
    return x, w


def jacobi_polynomial(x: npt.NDArray, n: int, alpha: float, beta: float) -> npt.NDArray[np.float64]:
    """
    Evaluates first `n+1` Jacobi polynomials at points `x` with parameters `alpha` and `beta`.
    Reflects L2, slide 12.

    Note that if multiple orders are desired, it would be much more efficient to
    simply return the full array P, but this function is designed to match the
    mathematical notation as closely as possible.

    Arguments:
        x: Points at which to evaluate the polynomials, shape (m,)
        n: Order of polynomial to compute (must be positive)
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$

    Returns: Array of shape (m,)
    """
    assert n >= 0, "n must be non-negative"

    P = np.empty((len(x), n+1))

    P[:, 0] = 1.0
    if n == 0:
        return P[:, n]

    P[:, 1] = 1 / 2 * (alpha - beta + (alpha + beta + 2) * x)
    if n == 1:
        return P[:, n]

    for k in range(1, n):
        a_nm1_n = (
            2
            * (k + alpha)
            * (k + beta)
            / ((2 * k + alpha + beta + 1) * (2 * k + alpha + beta))
        )
        a_n_n = (alpha**2 - beta**2) / (
            (2 * k + alpha + beta + 2) * (2 * k + alpha + beta)
        )
        a_np1_n = (
            2
            * (k + 1)
            * (k + alpha + beta + 1)
            / ((2 * k + alpha + beta + 2) * (2 * k + alpha + beta + 1))
        )

        P[:, k + 1] = ((a_n_n + x) * P[:, k] - a_nm1_n * P[:, k - 1]) / a_np1_n

    return P[:, n]


def jacobi_normalisation_const(
    n: int, alpha: float, beta: float
) -> int | npt.NDArray:
    """
    Computes the normalisation constant for Jacobi polynomials.
    Reflects $γ_n^(α,β)$ from L2, slide 11.

    Arguments:
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$
        n: int, $n ≥ 0$

    Returns: normalisation constant $γ_n$ for Jacobi polynomials $P_n$
    """
    return (
        2 ** (alpha + beta + 1)
        * (gamma(n + alpha + 1) * gamma(n + beta + 1))
        / (factorial(n) * (2 * n + alpha + beta + 1) * gamma(n + alpha + beta + 1))
    )


def jacobi_polynomial_normalised(
    x: npt.NDArray, n: int, alpha: float, beta: float
) -> npt.NDArray:
    """
    Convenience function to get normalized Jacobi polynomials.
    """
    P = jacobi_polynomial(x, n, alpha, beta)
    norm_const = jacobi_normalisation_const(n, alpha, beta)
    return P / np.sqrt(norm_const)

def grad_jacobi_polynomial(x: npt.NDArray, n: int, alpha: float, beta: float) -> npt.NDArray:
    """
    Computes the gradient of the first `n+1` Jacobi polynomials at nodes `x`.
    Reflects L2, slide 15.

    Arguments:
        x: Points at which to evaluate the gradients, shape (m,)
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$
        n: Highest order polynomial to compute (must be positive)
    Returns: Array of shape (m,) where each column corresponds to the gradient of a Jacobi polynomial
    """
    if n == 0:
        return np.zeros_like(x)

    coeff = 1 / 2 * (alpha + beta + n + 1)
    p_i = jacobi_polynomial(x, n=n - 1, alpha=alpha + 1, beta=beta + 1)

    return coeff * p_i

def grad_jacobi_polynomial_normalised(x, n, k, alpha, beta):
    """
    Computes the `k`th gradient of the Jacobi polynomials of order `n` at nodes `x`.
    Reflects L2, slide 21.

    Arguments:
        x: Points at which to evaluate the gradients, shape (m,)
        n: Order of polynomial to compute (must be positive)
        k: Order of derivative to compute (must be non-negative and less than or equal to n)
        alpha: Jacobi parameter, $α > -1$
        beta: Jacobi parameter, $β > -1$
    """
    if n == 0:
        return np.zeros_like(x)

    coeff = (
        1
        * gamma(alpha + beta + n + 1 + k)
        / (2**k * gamma(alpha + beta + n + 1))
        * np.sqrt(
            jacobi_normalisation_const(n-k, alpha + k, beta + k)
            / jacobi_normalisation_const(n, alpha, beta)
        )
    )
    p = jacobi_polynomial_normalised(x, n - k, alpha + k, beta + k)

    return coeff * p


class BasisFunction(Protocol):
    def __call__(self, x: npt.NDArray[np.float64], n: int) -> npt.NDArray[np.float64]:
        ...

def construct_vandermonde(
    x: npt.NDArray[np.float64], N: int, basis_function: BasisFunction
) -> npt.NDArray[np.float64]:
    """
    Construct the Vandermonde matrix V such that V[i,j] = phi_j(x_i),
    where phi_j is the j'th basis function evaluated at point x_i.
    """
    V = np.zeros((len(x), N), dtype=float)
    for j in range(N):
        V[:, j] = basis_function(x, j)
    return V


