Questions for Allan:
- Why subtract 0.5 flux? Chat says it comes from L operator
- We only converge to 10^-2
- How should we impose boundary conditions? Strongly or weakly? We get boundary integral, so weakly, but then it does not match lectures or FEM 5.2?


- Governing equation (is there a document or paper we should use as reference?)
- Could he elaborate a bit on first point, hints on rewrite of mathematical model to get skew symmetry
- Hints/tools to do the mesh for a circular domain
  - Arises from the fact that rotation matrices are skew symmetric?
- What exactly is the difference between continuous and discontinuous Galerkin?
  - Intuition: how can we get $K N^p$ scaling if we impose continuous BC's between elements? Does this not couple the points of elements together?
- Should we use continuous or discontinuous method?

Hints:
- Upwind boundary conditions
- Create grid for each element
- Global indices for nodes in grid

For week 12:
- Connectivity grid
- Global index

