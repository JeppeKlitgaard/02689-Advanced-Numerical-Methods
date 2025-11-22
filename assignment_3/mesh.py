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
        node_to_idx=node_to_idx,
    )
