import numpy as np
from scipy.special import eval_jacobi, factorial, gamma
from dataclasses import dataclass
from operators import Vandermonde2D_orthonormal, GradVandermonde2D_orthonormal, Dmatrices2D
import numpy.typing as npt
from mesh import Mesh

# quesiton
# 1) should eval_jacobi (NOT orthonormal) be changed to JacobiP (orthonormal)?

@dataclass(frozen=True, slots=True)
class DiscretizationElement:
    N: int  # Polynomial order N
    Mp: int # Number of nodes per element
    r: npt.NDArray[np.float64]
    s: npt.NDArray[np.float64]
    LocalReorder: npt.NDArray[np.int64]
    fmask_list: list[npt.NDArray[np.int64]]

    V: npt.NDArray[np.float64]
    Vr: npt.NDArray[np.float64]
    Vs: npt.NDArray[np.float64]
    Dr: npt.NDArray[np.float64]
    Ds: npt.NDArray[np.float64]
    M_canonical: npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class DiscretizationMesh:
    C: npt.NDArray[np.float64]
    gidx: npt.NDArray[np.int64]  

    x_full: npt.NDArray[np.float64]
    y_full: npt.NDArray[np.float64]
    x_global: npt.NDArray[np.float64]
    y_global: npt.NDArray[np.float64]

    BC_nodes: npt.NDArray[np.int64]


def create_discretization_element(N: int, LocalReorder = None) -> DiscretizationElement:

    x_equilateral, y_equilateral = Nodes2D(N)
    r_unordered, s_unordered = equilateral2rs(x_equilateral, y_equilateral)
    
    Mp = (N + 1) * (N + 2) // 2
    fmask1, fmask2, fmask3 = face_masks_triangle(r_unordered, s_unordered, tol=1e-10)
    if LocalReorder is None:
        LocalReorder = construct_LocalReorder(Mp, N, fmask1, fmask2, fmask3)

    r, s = r_unordered[LocalReorder], s_unordered[LocalReorder]

    # here the r, s are already reordered
    # V = Vandermonde2D_orthonormal(N, r_unordered, s_unordered, LocalReorder=LocalReorder)
    # Vr, Vs = GradVandermonde2D_orthonormal(N, r_unordered, s_unordered, LocalReorder=LocalReorder)
    V = Vandermonde2D_orthonormal(N, r, s)
    Vr, Vs = GradVandermonde2D_orthonormal(N, r, s)

    Dr, Ds = Dmatrices2D(V=V, Vr=Vr, Vs=Vs)
    M_canonical = np.linalg.inv(V @ V.T)

    return DiscretizationElement(
        N=N,
        Mp=Mp,
        r=r,
        s=s,
        LocalReorder=LocalReorder,
        V=V,
        Vr=Vr,
        Vs=Vs,
        Dr=Dr,
        Ds=Ds,
        M_canonical=M_canonical,
        fmask_list = [fmask1, fmask2, fmask3]
    )

def create_discretization_mesh(mesh: Mesh, discretization_element: DiscretizationElement) -> DiscretizationMesh:

    r, s = discretization_element.r, discretization_element.s
    x_full, y_full = construct_xy_full(mesh, r, s)
    
    C, gidx = build_global_map(discretization_element.N, mesh.EtoV, mesh.EtoE, mesh.EtoF)
    x_global, y_global = build_global_coords(C.T, x_full, y_full)

    fmask_list = discretization_element.fmask_list
    LocalReorder = discretization_element.LocalReorder
    BC_nodes, _ = boundary_nodes_from_connectivity(mesh.EtoE, C, fmask_list, LocalReorder)

    return DiscretizationMesh(
        C=C,
        gidx=gidx,
        x_full=x_full,
        y_full=y_full,
        x_global=x_global,
        y_global=y_global,
        BC_nodes=BC_nodes
    )


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


def Psi(r, s, v1, v2, v3):
    x = -0.5 * (r + s) * v1[0] + 0.5 * (1 + r) * v2[0] + 0.5 * (1 + s) * v3[0]
    y = -0.5 * (r + s) * v1[1] + 0.5 * (1 + r) * v2[1] + 0.5 * (1 + s) * v3[1]
    return x, y


def equilateral2rs(x, y):
    L1 = (np.sqrt(3)*y + 1) / 3
    L2 = (-3*x - np.sqrt(3)*y + 2) / 6
    L3 = (3*x - np.sqrt(3)*y + 2) / 6
    r = -L2 + L3 - L1
    s = -L2 - L3 + L1
    return r, s


# Chat fully
def warpfactor(N, rout, alpha=0, beta=0):
    """
    Warp factor (Warp–Blend) for order N evaluated at rout (array in [-1,1]).
    """
    rout = np.asarray(rout, dtype=float)          # size Nr

    # 1) LGL (Jacobi–Gauss–Lobatto) nodes
    r_GL = JacobiGL(alpha=alpha, beta=beta, N=N)  # size N+1

    # 2) Equidistant nodes
    req = np.linspace(-1.0, 1.0, N+1)

    # 3) Vandermonde at equidistant nodes: Veq[j,i] = P_i(req[j])
    Veq = np.column_stack([eval_jacobi(i, alpha, beta, req) for i in range(N+1)])  # (N+1)x(N+1)

    # 4) Evaluate basis at rout (Pmat[j,i] = P_i(rout[j]))
    Pmat = np.column_stack([eval_jacobi(i, alpha, beta, rout) for i in range(N+1)])  # (Nr)x(N+1)

    # 5) Lagrange basis values at rout on req nodes:
    # Veq^T Lmat^T = Pmat^T  => Lmat^T = (Veq^T)^{-1} Pmat^T
    Lmat_T = np.linalg.solve(Veq.T, Pmat.T)       # (N+1)xNr
    # Each column j: ℓ_i(rout_j)

    # 6) Displacement d_i
    disp = r_GL - req                             # (N+1,)

    # 7) Interpolate displacement
    warp = (Lmat_T.T @ disp)                      # (Nr,)

    # 8) Blend (vanish at endpoints)
    zerof = (np.abs(rout) < 1.0 - 1.0e-10)
    sf = 1.0 - (zerof * rout)**2                  # (Nr,)
    warp = warp / sf + warp * (~zerof) * (-1.0)   # when |r|≈1 force to 0

    return warp


# Chat fully
def Nodes2D(N):
    """
    Compute (x,y) warp–blend nodes in an equilateral triangle for polynomial order N.
    Returns arrays x,y of length (N+1)(N+2)/2.
    """
    if N == 0:
        # Triangle vertices (equilateral of side 2 centered)
        return np.array([0.0]), np.array([0.0])
    # Optimized alpha (same list as in original code; MATLAB 1-based → Python 0-based)
    alpopt = [0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999,
              1.2832, 1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258]
    alpha = alpopt[N] if N < 15 else 5/3

    Np = (N+1)*(N+2)//2
    L1 = np.zeros(Np); L2 = np.zeros(Np); L3 = np.zeros(Np)
    sk = 0
    for n in range(1, N+2):
        for m in range(1, N+3 - n):
            L1[sk] = (n-1)/N
            L3[sk] = (m-1)/N
            sk += 1
    L2 = 1.0 - L1 - L3

    # Initial equidistributed coordinates in equilateral triangle
    x = -L2 + L3
    y = (-L2 - L3 + 2*L1)/np.sqrt(3.0)

    # Blending functions
    blend1 = 4.0 * L2 * L3
    blend2 = 4.0 * L1 * L3
    blend3 = 4.0 * L1 * L2

    # Edge warp amounts
    warpf1 = warpfactor(N, L3 - L2)
    warpf2 = warpfactor(N, L1 - L3)
    warpf3 = warpfactor(N, L2 - L1)

    # Combine (edge-wise scaling factors)
    warp1 = blend1 * warpf1 * (1.0 + (alpha * L1)**2)
    warp2 = blend2 * warpf2 * (1.0 + (alpha * L2)**2)
    warp3 = blend3 * warpf3 * (1.0 + (alpha * L3)**2)

    # Accumulate deformations (rotations 0, 120°, 240°)
    x = x + 1.0*warp1 + np.cos(2*np.pi/3)*warp2 + np.cos(4*np.pi/3)*warp3
    y = y + 0.0*warp1 + np.sin(2*np.pi/3)*warp2 + np.sin(4*np.pi/3)*warp3

    return x, y


# Allan FEM book: Algo 14
def build_global_map(P, EtoV, EtoE, EToF):
    EtoV = np.asarray(EtoV, dtype=int)
    EtoE = np.asarray(EtoE, dtype=int)
    EToF = np.asarray(EToF, dtype=int)
    K = EtoV.shape[0]
    Nfaces = 3
    MP = (P + 1) * (P + 2) // 2
    Mpf = P + 1
    face_interior_count = Mpf - 2          # P-1
    interior_count = (P - 1) * (P - 2) // 2
    Nv = EtoV.max() + 1                    # zero-based vertex count
    C = np.zeros((K, MP), dtype=int)
    gidx = Nv
    for n in range(K):
        C[n, :Nfaces] = EtoV[n, :Nfaces]
        # faces
        for i in range(Nfaces):
            start = Nfaces + i * face_interior_count
            end = start + face_interior_count
            if face_interior_count > 0:
                if EtoE[n, i] >= n:  # new or boundary
                    C[n, start:end] = np.arange(gidx, gidx + face_interior_count)
                    gidx += face_interior_count
                else:
                    k_neighbor = EtoE[n, i]
                    f_neighbor = EToF[n, i]
                    n_start = Nfaces + f_neighbor * face_interior_count
                    n_end = n_start + face_interior_count
                    C[n, start:end] = C[k_neighbor, n_start:n_end][::-1]
        # interior
        if interior_count > 0:
            int_start = Nfaces + Nfaces * face_interior_count
            int_end = MP
            C[n, int_start:int_end] = np.arange(gidx, gidx + interior_count)
            gidx += interior_count
    return C, gidx


def build_global_coords(C, x_elem, y_elem):
    # x_elem, y_elem: list/array of shape (K, Mp) with element-local coords
    K, Mp = C.shape
    gidx = int(C.max()) + 1
    xg = np.full(gidx, np.nan)
    yg = np.full(gidx, np.nan)
    for k in range(K):
        g = C[k]               # global ids for element k
        xk = x_elem[k]         # (Mp,)
        yk = y_elem[k]
        # scatter with consistency check
        for i, gi in enumerate(g):
            if np.isnan(xg[gi]):
                xg[gi] = xk[i]; yg[gi] = yk[i]
            else:
                # shared nodes should match
                assert np.allclose([xg[gi], yg[gi]], [xk[i], yk[i]])
    return xg, yg


def apply_local_reorder_to_faces(faces_local, LocalReorder=None):
    """
    Map face masks (built in canonical local ordering) to a reordered local DOF ordering.
    faces_local: list of 3 arrays (local node indices for each face)
    LocalReorder: 1D permutation of length Mp (new_order = old[LocalReorder]).
    """
    if LocalReorder is None:
        return faces_local
    LocalReorder = np.asarray(LocalReorder, dtype=int)
    inv = np.empty_like(LocalReorder)
    inv[LocalReorder] = np.arange(LocalReorder.size)
    return [inv[np.asarray(f, dtype=int)] for f in faces_local]

def boundary_nodes_from_connectivity(EToE, C, faces_local, LocalReorder=None):
    """
    Collect ONLY domain-boundary global DOF indices (no interior edges).

    EToE: (K,3) element-to-element connectivity
          boundary faces are marked as self (e) or sentinel (0 or -1), depending on convention.
    C:    (K,Mp) local-to-global map (0-based)
    faces_local: list of 3 arrays of local node indices per face, matching EToE's local face order
                 e.g. [(1,2),(2,3),(3,1)] masks built on your reference nodes.
    LocalReorder: optional 1D permutation of local DOF ordering used in C and element operators.

    Returns:
      bd_nodes: unique sorted global DOF indices on the outer boundary
      bd_faces: list of arrays with global DOFs per boundary face occurrence (one per boundary face)
    """
    EToE = np.asarray(EToE, dtype=int)
    C    = np.asarray(C, dtype=int)
    K, Nfaces = EToE.shape
    assert Nfaces == 3

    # map face masks if local DOFs were reordered
    faces_local = apply_local_reorder_to_faces(faces_local, LocalReorder)

    # Determine boundary faces across common conventions
    eidx = np.arange(K)[:, None]
    if EToE.min() >= 1:       # 1-based elements, boundary often as self=e+1 or 0
        is_self = (EToE == (eidx + 1))
        is_bdy = is_self | (EToE == 0)
    else:                     # 0-based elements, boundary often as self=e or -1
        is_self = (EToE == eidx)
        is_bdy = is_self | (EToE == -1)

    bd_faces = []
    for e in range(K):
        for f in range(3):
            if is_bdy[e, f]:
                bd_faces.append(C[e, faces_local[f]])

    bd_nodes = np.unique(np.concatenate(bd_faces)) if bd_faces else np.array([], dtype=int)
    return bd_nodes, bd_faces


# def prepare_general_operators(N, LocalReorder):

#     x_equilateral, y_equilateral = Nodes2D(N)
#     r, s = equilateral2rs(x_equilateral, y_equilateral)
#     r_reordered = r[LocalReorder]
#     s_reordered = s[LocalReorder]
#     V = Vandermonde2D_orthonormal(N, r, s, LocalReorder=LocalReorder)
#     Vr, Vs = GradVandermonde2D_orthonormal(N, r, s, LocalReorder=LocalReorder)
#     Dr, Ds = Dmatrices2D(V=V, Vr=Vr, Vs=Vs)
#     M_inv_canonical = V @ V.T
#     M_canonical = np.linalg.inv(M_inv_canonical)
    
#     return r_reordered, s_reordered, Dr, Ds, M_canonical


def face_masks_triangle(r, s, tol=1e-10):
    fmask1 = np.where(np.abs(s + 1.0) < tol)[0]          # s = -1
    fmask2 = np.where(np.abs(r + s) < tol)[0]            # r + s = 0
    fmask3 = np.where(np.abs(r + 1.0) < tol)[0]          # r = -1
    return fmask1, fmask2, fmask3


def construct_LocalReorder(Mp, N, fmask1, fmask2, fmask3):
    Fmask = np.concatenate([fmask1, fmask2, fmask3])
    f_interior = np.setdiff1d(np.arange(0, Mp-1), Fmask)
    LocalReorder = np.concatenate([
        np.array([0, N, Mp -1]),
        fmask1[1:N],            # fid1(2:Mpf-1)
        fmask2[1:N],            # fid2(2:Mpf-1)
        fmask3[N-1:0:-1],       # fid3(Mpf-1:-1:2)
        f_interior
    ])
    return LocalReorder


def construct_xy_full(mesh: Mesh, r, s):
    R, S = r[:, None], s[:, None]  # (Np,1)
    v1, v2, v3 = mesh.EtoV[:, 0], mesh.EtoV[:, 1], mesh.EtoV[:, 2]
    V_x, V_y = mesh.V_x, mesh.V_y
    x_full = 0.5 * (-(R + S) * V_x[v1] + (1.0 + R) * V_x[v2] + (1.0 + S) * V_x[v3])
    y_full = 0.5 * (-(R + S) * V_y[v1] + (1.0 + R) * V_y[v2] + (1.0 + S) * V_y[v3])
    return x_full, y_full