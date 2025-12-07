import scipy.sparse as sp
import numpy as np
import scipy.sparse as sp
import numpy.typing as npt
from tqdm import tqdm
from solver import get_n_element_operators

class SolverHeatEquation2D:
    def __init__(
        self,
        mesh,
        discretization_element,
        discretization_mesh,
        advection_velocity_field,
        g_bc,
        g_bc_dt,
        f_rhs,
        alpha
    ):
        self.mesh = mesh
        self.discretization_element = discretization_element
        self.discretization_mesh = discretization_mesh
        self.advection_velocity_field = advection_velocity_field
        self.g_bc = g_bc  # Boundary Condition u(x,y,t)
        self.g_bc_dt = g_bc_dt  # Time derivative of BC du/dt(x,y,t)
        self.f_rhs = f_rhs  # Source function f(x,y,t)
        self.alpha = alpha
        self.L_global, self.M_global = self.construct_global_assembly()
        self.M_sys, self.L_sys = self.apply_matrix_bc()
        self.M_sys_solver = sp.linalg.factorized(self.M_sys.tocsc())

    def construct_L_n_element(self, Dx, Dy, M, x_nodes, y_nodes):
        return self.alpha * (Dx.T @ M @ Dx + Dy.T @ M @ Dy)

    def construct_global_assembly(self) -> tuple[npt.NDArray, npt.NDArray]:
        """
        Implementation of Global Assembly based on FEM Book Algorithms 15 & 16.
        Modified for Advection (Non-Symmetric L).
        """
        # print("Assembling global matrices...")

        gidx = self.discretization_mesh.gidx

        # Use sparse matrix construction for efficiency
        # We collect triplets (row, col, value)
        L_rows, L_cols, L_vals = [], [], []
        M_rows, M_cols, M_vals = [], [], []

        C = self.discretization_mesh.C
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global

        # for n_element in tqdm(range(self.mesh.num_elements)):
        for n_element in range(self.mesh.num_elements):
            # 1. Get Element Geometry
            # (Assuming get_n_element_operators is available as in your snippet)
            # You might need to pass physical nodes to your operator function if it calculates metrics
            node_indices = C[n_element, :]
            x_nodes = x_global[node_indices]
            y_nodes = y_global[node_indices]

            # Reconstruct vertices for the metric terms if needed by your helper
            # (Assuming your existing helper works as before)
            x_vertex = self.mesh.V_x[self.mesh.EtoV[n_element, :]]
            y_vertex = self.mesh.V_y[self.mesh.EtoV[n_element, :]]
            v1 = np.array([x_vertex[0], y_vertex[0]])
            v2 = np.array([x_vertex[1], y_vertex[1]])
            v3 = np.array([x_vertex[2], y_vertex[2]])

            # Get local operators (Reference Dx, Dy and Metric-scaled M)
            Dx, Dy, M_local, _, _ = get_n_element_operators(
                discretization_element=self.discretization_element,
                v1_n=v1,
                v2_n=v2,
                v3_n=v3,
            )

            # Construct Local Advection Matrix L (Algorithm 15)
            L_local = self.construct_L_n_element(Dx, Dy, M_local, x_nodes, y_nodes)

            # 3. Add to lists for Sparse Assembly
            # Create meshgrid of indices for this element
            r_idx, c_idx = np.meshgrid(node_indices, node_indices, indexing="ij")

            L_rows.extend(r_idx.flatten())
            L_cols.extend(c_idx.flatten())
            L_vals.extend(L_local.flatten())

            M_rows.extend(r_idx.flatten())
            M_cols.extend(c_idx.flatten())
            M_vals.extend(M_local.flatten())

        # Create Sparse Matrices
        L_global = sp.coo_matrix((L_vals, (L_rows, L_cols)), shape=(gidx, gidx)).tocsr()
        M_global = sp.coo_matrix((M_vals, (M_rows, M_cols)), shape=(gidx, gidx)).tocsr()

        return L_global, M_global


    def apply_matrix_bc(self):
        M_sys = self.M_global.copy()
        L_sys = self.L_global.copy()
        bc_nodes = self.discretization_mesh.BC_nodes

        M_sys[:, bc_nodes] = 0.0  # Zero out columns
        M_sys[bc_nodes, :] = 0.0  # Zero out rows
        M_sys[bc_nodes, bc_nodes] = 1.0  # Set diagonal to 1.0

        # Operator Matrix
        L_sys[:, bc_nodes] = 0.0  # Zero out columns
        L_sys[bc_nodes, :] = 0.0  # Zero out rows

        return M_sys, L_sys

    def get_rhs_IVP(self, t, u):
        """
        Computes du/dt = M^-1 * (-L*u + f + BC_terms)
        """
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global
        bc_nodes = self.discretization_mesh.BC_nodes

        u_full = u.copy()
        u_full[bc_nodes] = self.g_bc(x_global[bc_nodes], y_global[bc_nodes], t)

        # Construct RHS
        # Term 1: - L * u
        rhs = -self.L_sys @ u_full
        # Term 2: + f (source term)
        f_vec = self.f_rhs(x_global, y_global, t)
        rhs += self.M_sys @ f_vec

        rhs[bc_nodes] = self.g_bc_dt(x_global[bc_nodes], y_global[bc_nodes], t)

        du_dt = self.M_sys_solver(rhs)

        return du_dt