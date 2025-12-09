from discretization import DiscretizationElement, Psi
from operators import GeometricFactors2D, Dmatrices2D_xy
import numpy as np
import scipy.sparse as sp
import numpy.typing as npt
from tqdm import tqdm
from utils import jacobi_gauss_lobatto

# this could be computed in DiscretizationMesh class
def get_n_element_operators(discretization_element: DiscretizationElement, v1_n, v2_n, v3_n):

    x_n, y_n = Psi(discretization_element.r, discretization_element.s, v1_n, v2_n, v3_n)
    rx, sx, ry, sy, J = GeometricFactors2D(x_n, y_n, discretization_element.Dr, discretization_element.Ds)
    M = J[:, None] * discretization_element.M_canonical
    Dx, Dy = Dmatrices2D_xy(Dr=discretization_element.Dr, Ds=discretization_element.Ds, rx=rx, sx=sx, ry=ry, sy=sy)

    return Dx, Dy, M, x_n, y_n


class AdvectionIVPSolver2D:
    def __init__(
        self,
        mesh,
        discretization_element,
        discretization_mesh,
        advection_velocity_field,
        g_bc,
        g_bc_dt,
        f_rhs,
        with_upwind_scheme: bool = True,
    ):
        self.mesh = mesh
        self.discretization_element = discretization_element
        self.discretization_mesh = discretization_mesh
        self.with_upwind_scheme = with_upwind_scheme

        # Advection velocity field, must be stationary
        self.advection_velocity_field = advection_velocity_field

        # Resolve boundary nodes to impose boundary conditions at while respecting upwind scheme
        self.resolved_boundary_nodes = self.resolve_boundary_nodes()

        # These are used for the boundary conditions
        self.g_bc = g_bc  # Boundary Condition u(x,y,t)
        self.g_bc_dt = g_bc_dt  # Time derivative of BC du/dt(x,y,t)

        # This is the source term/RHS without
        self.f_rhs = f_rhs  # Source function f(x,y,t)

        # Assemble constant matrices immediately (since velocity field is stationary)
        self.L_global, self.M_global = self.construct_global_assembly()

        # Apply BCs to get system matrices
        self.M_sys, self.L_sys = self.apply_matrix_bc()

        # Pre-factorize M for efficiency
        self.M_sys_solver = sp.linalg.factorized(self.M_sys.tocsc())


    def resolve_boundary_nodes(self) -> npt.NDArray:
        """
        Identify boundary nodes where the advection velocity points outward.
        These nodes will have their values imposed directly from the boundary condition.
        """
        if not self.with_upwind_scheme:
            return self.discretization_mesh.BC_nodes

        bc_nodes = self.discretization_mesh.BC_nodes
        bc_normals = self.discretization_mesh.BC_normals

        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global

        a_x, a_y = self.advection_velocity_field(x_global[bc_nodes], y_global[bc_nodes])

        a_n = a_x * bc_normals[:, 0] + a_y * bc_normals[:, 1]

        return bc_nodes[a_n < 0]

    def construct_L_n_element(self, Dx, Dy, M, x_nodes, y_nodes):
        ax, ay = self.advection_velocity_field(x_nodes, y_nodes)

        # Conservative form
        # Strong form advection operator on element: M @ (ax*Dx + ay*Dy)
        # Corresponds to (v, a.grad(u))
        # return M @ (np.diag(ax) @ Dx + np.diag(ay) @ Dy)

        # Split form
        MAx = M @ np.diag(ax) @ Dx
        MAy = M @ np.diag(ay) @ Dy
        return 1/2 * ((MAx - MAx.T) + (MAy - MAy.T))

    def construct_global_assembly(self) -> tuple[npt.NDArray, npt.NDArray]:
        """
        Implementation of Global Assembly based on FEM Book Algorithms 15 & 16.
        Modified for Advection (Non-Symmetric L).
        """
        print("Assembling global matrices...")

        gidx = self.discretization_mesh.gidx

        # Use sparse matrix construction for efficiency
        # We collect triplets (row, col, value)
        L_rows, L_cols, L_vals = [], [], []
        M_rows, M_cols, M_vals = [], [], []

        C = self.discretization_mesh.C
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global

        for n_element in tqdm(range(self.mesh.num_elements)):
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

    # def apply_matrix_bc(self):
    #     """
    #     Modifies M and L to enforce BCs on the system level (Algorithm 17).
    #     For boundary node i:
    #       Row i of L_sys becomes 0.
    #       Row i of M_sys becomes unit row (diagonal=1, others=0).
    #     Equation becomes: 1 * du_i/dt = RHS_i
    #     """
    #     M_sys = self.M_global.copy()
    #     L_sys = self.L_global.copy()
    #     bc_nodes = self.resolved_boundary_nodes

    #     # Zero out rows for boundary nodes
    #     # Note: We do NOT zero out columns. The interior nodes still need
    #     # to "feel" the mass/advection from the boundary nodes.
    #     M_sys[bc_nodes, :] = 0.0
    #     L_sys[bc_nodes, :] = 0.0

    #     # Set diagonal of Mass matrix to 1.0 for boundary nodes
    #     M_sys[bc_nodes, bc_nodes] = 1.0

    #     return M_sys, L_sys

    def apply_matrix_bc(self):
        """
        Modifies M and L to enforce BCs on the system level (Algorithm 17).
        For boundary node i:
          Row i of L_sys becomes 0.
          Row i of M_sys becomes unit row (diagonal=1, others=0).
        Equation becomes: 1 * du_i/dt = RHS_i
        """

        # IMPOSE WEAKLY INSTEAD!
        return self.M_global, self.L_global

    def compute_boundary_flux_term(self, u_current, t):
        """
        Computes the boundary integral vector using the Upwind Flux from Eq (16).
        """
        flux_vector = np.zeros(self.discretization_mesh.gidx)
        _, w_1d = jacobi_gauss_lobatto(self.discretization_element.N, 0, 0)

        boundary_faces = self.discretization_mesh.BC_faces
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global

        for face_nodes in boundary_faces:
            x_f = x_global[face_nodes]
            y_f = y_global[face_nodes]

            # Calculate Edge Length
            dx = x_f[-1] - x_f[0]
            dy = y_f[-1] - y_f[0]
            length = np.sqrt(dx**2 + dy**2)
            J_s = length / 2.0

            # Normal vector components
            nx = self.discretization_mesh.BC_normals[face_nodes[0], 0]
            ny = self.discretization_mesh.BC_normals[face_nodes[0], 1]

            # Velocity at face nodes
            ax, ay = self.advection_velocity_field(x_f, y_f)
            a_n = ax * nx + ay * ny

            # State values
            u_minus = u_current[face_nodes]       # Internal solution
            u_bc = self.g_bc(x_f, y_f, t)         # Boundary condition

            # If a_n > 0 (Outflow): f* = a_n * u_minus
            # If a_n < 0 (Inflow):  f* = a_n * u_bc
            f_n_star = np.where(a_n > 0,
                                a_n * u_minus,
                                a_n * u_bc)

            net_flux = f_n_star - (a_n * u_minus)

            # Integrate over the face
            weighted_flux = net_flux * w_1d * J_s
            flux_vector[face_nodes] += weighted_flux

        return flux_vector

    def get_rhs_IVP(self, t, u):
        """
        Computes du/dt = M^-1 * (-L*u + f + BC_terms)
        """
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global

        # Volume contribution
        rhs = -self.L_global @ u

        # Source term
        f_vals = self.f_rhs(x_global, y_global, t)
        rhs += self.M_global @ f_vals

        # Boundary flux
        flux = self.compute_boundary_flux_term(u, t)
        rhs -= flux

        # Solve
        du_dt = self.M_sys_solver(rhs)

        return du_dt