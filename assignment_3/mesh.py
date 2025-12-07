"""
Contains basic meshing functionality.
"""

import pickle
from dataclasses import asdict, dataclass
from typing import Self

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.sparse import coo_matrix, eye


@dataclass(frozen=True, slots=True)
class Mesh:
    num_elements: int  # K in Allan's notation
    num_nodes: int  # Nv in Allan's notation
    num_faces: int  # Nfaces in Allan's notation
    V_x: npt.NDArray[np.float64]  # Node x-coordinates, VX in Allan's notation
    V_y: npt.NDArray[np.float64]  # Node y-coordinates, VY in Allan's notation

    # Element to vertex connectivity
    # Elements this array are node_tags, NOT INDICES!
    EtoV_tags: npt.NDArray[
        np.int64
    ]  # Element to vertex connectivity, EToV in Allan's notation
    # Elements in this array are indices (zero-based)
    EtoV: npt.NDArray[
        np.int64
    ]  # Element to vertex connectivity, EToV in Allan's notation
    EtoE: npt.NDArray[
        np.int64
    ]  # comment here
    EtoF: npt.NDArray[
        np.int64
    ]  # comment here

    # Node-index mapping
    node_to_idx: dict[int, int]  # Maps node tags to zero-based indices

    def save(self, filename: str) -> None:
        """
        Save the mesh to a file using pickle.

        File should be named something like `XXX.mesh.pkl`

        Parameters:
            filename: The name of the file to save the mesh to.
        """
        with open(filename, "wb") as f:
            pickle.dump(asdict(self), f)

    @classmethod
    def load(cls, filename: str) -> Self:
        """
        Load a mesh from a file.

        Parameters:
            filename: The name of the file to load the mesh from.
        Returns:
            Mesh: The loaded mesh object.
        """
        with open(filename, "rb") as f:
            return cls(**pickle.load(f))

    def plot(self, ax=None) -> None:
        """
        Plots the mesh using matplotlib.

        Can optionally provide an axis to plot on.

        Parameters:
            ax: Matplotlib axis to plot on. If None, uses current axis.
        """
        if ax is None:
            ax = plt.gca()
            ax.set_aspect("equal")
        ax.triplot(self.V_x, self.V_y, triangles=self.EtoV, color="blue")
        ax.set_aspect("equal")


def create_2d_circle(
    radius: float, origin: tuple[float, float], mesh_size: float
) -> Mesh:
    """
    Create a 2D circle and mesh it using triangles via Gmsh.

    Parameters:
        radius: Radius of the circle.
        origin: (x, y) coordinates of the circle's center.
        mesh_size: Maximum size of the mesh elements.
            Larger values do not necessarily yield fewer elements
            past a certain point

    Returns:
        Mesh: The generated mesh object.
    """
    origin_3d = (origin[0], origin[1], 0.0)

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

    # Create circle geometry
    gmsh.model.add("circle")
    tag_circle = gmsh.model.occ.addCircle(*origin_3d, radius)
    tag_loop = gmsh.model.occ.addCurveLoop([tag_circle])
    tag_surface = gmsh.model.occ.addPlaneSurface([tag_loop])

    # Synchronize CAD and Gmsh model
    gmsh.model.occ.synchronize()

    # Generate mesh
    gmsh.model.mesh.generate(dim=2)

    # Do some basic validation and sanity checks
    element_types, element_tags, node_tags_element = gmsh.model.mesh.getElements(
        dim=2, tag=tag_surface
    )
    assert len(element_types) == 1, "Unexpected number of element types"
    element_type = element_types[0]
    num_elements = len(element_tags[0])

    name, dim, order, num_nodes_per_element, local_node_coord, num_primary_nodes = (
        gmsh.model.mesh.getElementProperties(element_type)
    )
    assert dim == 2, "Element dimension mismatch"
    assert num_nodes_per_element == 3, "We want triangular mesh, something went wrong"

    # Construct node coordinate arrays
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    # node_tags, node_coords, _ = gmsh.model.mesh.getNodes(dim=2, tag=tag_surface)  # I don't understand why this is not correct, but I think it just excludes the boundary for reasons unknown.
    # node_coords is flattened: [x1, y1, z1, x2, y2, z2, ...]
    node_coords = np.array(node_coords).reshape(-1, 3)
    V_x = node_coords[:, 0]
    V_y = node_coords[:, 1]

    # Construct EtoV
    # nodeTags[0] contains the flattened list of node tags for all elements
    # We can reshape this using numpy for easier viewing
    # The list is [e1_n1, e1_n2, e1_n3, e2_n1, e2_n2, e2_n3, ...]
    node_to_idx = {tag: idx for idx, tag in enumerate(node_tags)}
    EtoV_tags = np.array(node_tags_element).reshape(num_elements, num_nodes_per_element)
    EtoV = np.vectorize(node_to_idx.__getitem__)(EtoV_tags)

    EtoE, EtoF = tiConnect2D(EtoV)

    num_nodes = len(node_tags)

    # Planar graph, so _obviouslyly_, number of vertices is simply:
    num_faces = 2 + num_elements - num_nodes

    gmsh.finalize()

    return Mesh(
        num_elements=num_elements,
        num_nodes=num_nodes,
        num_faces=num_faces,
        V_x=V_x,
        V_y=V_y,
        EtoV_tags=EtoV_tags,
        EtoV=EtoV,
        EtoE=EtoE,
        EtoF=EtoF,
        node_to_idx=node_to_idx,
    )


def create_2d_square(
    half_width: float, origin: tuple[float, float], mesh_size: float
) -> Mesh:
    """
    Create a 2D square and mesh it using triangles via Gmsh.

    Parameters:
        half_width: Half of the square's side length.
        origin: (x, y) coordinates of the square's center.
        mesh_size: Maximum size of the mesh elements.

    Returns:
        Mesh: The generated mesh object.
    """
    ox, oy = origin
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
    gmsh.model.add("square")

    # Square corner points in CCW order
    p1 = gmsh.model.occ.addPoint(ox - half_width, oy - half_width, 0)
    p2 = gmsh.model.occ.addPoint(ox + half_width, oy - half_width, 0)
    p3 = gmsh.model.occ.addPoint(ox + half_width, oy + half_width, 0)
    p4 = gmsh.model.occ.addPoint(ox - half_width, oy + half_width, 0)

    # Lines forming boundary
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)
    l4 = gmsh.model.occ.addLine(p4, p1)

    loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    surface = gmsh.model.occ.addPlaneSurface([loop])

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.generate(dim=2)

    # Get element data
    element_types, element_tags, node_tags_element = gmsh.model.mesh.getElements(
        dim=2, tag=surface
    )
    assert len(element_types) == 1, "Unexpected number of element types"
    element_type = element_types[0]
    num_elements = len(element_tags[0])

    name, dim, order, num_nodes_per_element, local_node_coord, num_primary_nodes = (
        gmsh.model.mesh.getElementProperties(element_type)
    )
    assert dim == 2, "Element dimension mismatch"
    assert num_nodes_per_element == 3, "We want triangular mesh, something went wrong"

    # Node coordinates
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = np.array(node_coords).reshape(-1, 3)
    V_x = node_coords[:, 0]
    V_y = node_coords[:, 1]

    # EtoV connectivity
    node_to_idx = {tag: idx for idx, tag in enumerate(node_tags)}
    EtoV_tags = np.array(node_tags_element).reshape(num_elements, num_nodes_per_element)
    EtoV = np.vectorize(node_to_idx.__getitem__)(EtoV_tags)

    EtoE, EtoF = tiConnect2D(EtoV)

    num_nodes = len(node_tags)

    # Euler characteristic formula for planar triangulation
    num_faces = 2 + num_elements - num_nodes

    gmsh.finalize()

    return Mesh(
        num_elements=num_elements,
        num_nodes=num_nodes,
        num_faces=num_faces,
        V_x=V_x,
        V_y=V_y,
        EtoV_tags=EtoV_tags,
        EtoV=EtoV,
        EtoE=EtoE,
        EtoF=EtoF,
        node_to_idx=node_to_idx,
    )



def tiConnect2D(EToV):

    EToV = np.asarray(EToV, dtype=int)
    K = EToV.shape[0]
    Nfaces = 3
    vn = np.array([[0,1],[1,2],[0,2]])  # local faces
    Nv = EToV.max()+1  # assume 0-based
    TotalFaces = K * Nfaces

    rows = []
    cols = []
    data = []
    sk = 0
    for k in range(K):
        for f in range(Nfaces):
            v0, v1 = EToV[k, vn[f]]
            rows.extend([sk, sk])
            cols.extend([v0, v1])
            data.extend([1, 1])
            sk += 1
    SpFToV = coo_matrix((data, (rows, cols)), shape=(TotalFaces, Nv)).tocsr()
    SpFToF = SpFToV @ SpFToV.T - 2 * eye(TotalFaces, format='csr')

    coo = SpFToF.tocoo()
    mask = coo.data == 2
    faces1 = coo.row[mask]
    faces2 = coo.col[mask]

    element1 = faces1 // Nfaces
    face1 = faces1 % Nfaces
    element2 = faces2 // Nfaces
    face2 = faces2 % Nfaces

    EToE = np.tile(np.arange(K)[:, None], (1, Nfaces))
    EToF = np.tile(np.arange(Nfaces)[None, :], (K, 1))
    for e1, f1, e2, f2 in zip(element1, face1, element2, face2):
        EToE[e1, f1] = e2
        EToF[e1, f1] = f2
    
    return EToE, EToF

def create_toy_mesh(V_x, V_y, EtoV):
    EtoE, EtoF = tiConnect2D(EtoV)
    num_elements = EtoV.shape[0]
    num_nodes = V_x.shape[0]
    num_faces = 2 + num_elements - num_nodes

    return Mesh(
        num_elements=num_elements,
        num_nodes=num_nodes,
        num_faces=num_faces,
        V_x=V_x,
        V_y=V_y,
        EtoV=EtoV,
        EtoE=EtoE,
        EtoF=EtoF,
        node_to_idx = dict(),
        EtoV_tags = None
    )