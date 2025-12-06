from discretization import DiscretizationElement, Psi
from operators import GeometricFactors2D, Dmatrices2D_xy

# this could be computed in DiscretizationMesh class
def get_n_element_operators(discretization_element: DiscretizationElement, v1_n, v2_n, v3_n):

    x_n, y_n = Psi(discretization_element.r, discretization_element.s, v1_n, v2_n, v3_n)
    rx, sx, ry, sy, J = GeometricFactors2D(x_n, y_n, discretization_element.Dr, discretization_element.Ds)
    M = J[:, None] * discretization_element.M_canonical
    Dx, Dy = Dmatrices2D_xy(Dr=discretization_element.Dr, Ds=discretization_element.Ds, rx=rx, sx=sx, ry=ry, sy=sy)

    return Dx, Dy, M, x_n, y_n