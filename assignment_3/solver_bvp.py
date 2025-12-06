import numpy as np
from abc import ABC, abstractmethod
from mesh import Mesh
from discretization import DiscretizationElement, DiscretizationMesh, Psi
from operators import Vandermonde2D_orthonormal, GradVandermonde2D_orthonormal, Dmatrices2D, Dmatrices2D_xy, GeometricFactors2D
from tqdm import tqdm

class BVP_solver(ABC):

    def __init__(self, mesh: Mesh, discretization_element: DiscretizationElement, discretization_mesh: DiscretizationMesh):
        self.mesh = mesh
        self.discretization_element = discretization_element
        self.discretization_mesh = discretization_mesh
        self.g_bc = None  # placeholder for boundary condition function
        self.f_rhs = None  # placeholder for RHS function
        self.L_N = None
        self.f_N = None
        self.u_N = None

    @abstractmethod
    def construct_L_n_element(self, Dx, Dy, M):
        pass 

    def construct_global_assembly(self):
        # construct L_N, f_N

        L_N = np.zeros((self.discretization_mesh.gidx, self.discretization_mesh.gidx))
        f_N = np.zeros(self.discretization_mesh.gidx)

        Mp = self.discretization_element.Mp
        C = self.discretization_mesh.C
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global

        # Allan FEM book: Algo 15, 16
        for n_element in tqdm(range(self.mesh.num_elements)):

            x_vertex_n = self.mesh.V_x[self.mesh.EtoV[n_element, :]]
            y_vertex_n = self.mesh.V_y[self.mesh.EtoV[n_element, :]]

            v1_n = np.array([x_vertex_n[0], y_vertex_n[0]])
            v2_n = np.array([x_vertex_n[1], y_vertex_n[1]])
            v3_n = np.array([x_vertex_n[2], y_vertex_n[2]])

            Dx, Dy, M, _, _ = get_n_element_operators(
                discretization_element=self.discretization_element,
                v1_n=v1_n, v2_n=v2_n, v3_n=v3_n
            )

            L_n_element = self.construct_L_n_element(Dx=Dx, Dy=Dy, M=M)

            for j in range(Mp):
                for i in range(Mp):
                    # exploit symmetry, not working
                    # if C[n_element, j] >= C[n_element, i]:
                    L_N[C[n_element, i], C[n_element, j]] += L_n_element[i, j]

            for j in range(Mp):
                jj = C[n_element, j]
                x_j = x_global[jj]
                y_j = y_global[jj]
                for i in range(Mp):
                    ii = C[n_element, i]
                    f_N[ii] += M[i, j] * self.f_rhs(x_j, y_j)

        self.L_N = L_N
        self.f_N = f_N


    def prescribe_boundary_conditions(self):

        BC_nodes = self.discretization_mesh.BC_nodes
        x_global = self.discretization_mesh.x_global
        y_global = self.discretization_mesh.y_global
        L_N = self.L_N.copy()
        f_N = self.f_N.copy()
        # Allan FEM book: Algo 17
        for i in BC_nodes:
            # correct: in FEM book the f function are KNOWN boundary conditions -> so u_exact
            u_i = self.g_bc(x_global[i], y_global[i])   # Dirichlet value
            f_N = f_N - L_N[:, i] * u_i      # adjust RHS
            L_N[i, :] = 0.0
            L_N[:, i] = 0.0
            L_N[i, i] = 1.0
            
        for i in BC_nodes:
            f_N[i] = self.g_bc(x_global[i], y_global[i])

        self.L_N = L_N
        self.f_N = f_N


    def solve(self):
        self.construct_global_assembly()
        self.prescribe_boundary_conditions()
        self.u_N = np.linalg.solve(self.L_N, self.f_N)


class BVP_Poisson_solver(BVP_solver):
    # test solver for the Poisson problem
    def __init__(self, mesh: Mesh, discretization_element: DiscretizationElement,
                 discretization_mesh: DiscretizationMesh, g_bc, f_rhs):
        super().__init__(mesh, discretization_element, discretization_mesh)
        self.g_bc = g_bc
        self.f_rhs = f_rhs

    def construct_L_n_element(self, Dx, Dy, M):
        return Dx.T @ M @ Dx + Dy.T @ M @ Dy

