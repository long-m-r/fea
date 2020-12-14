# Flux Envelope Analysis
A Python package and method for generating reduced-dimension solution spaces from large linear programs. It was designed for solving 2D and 3D solution spaces from large linear programs as utilized in constraint-based genome-scale metabolic modeling.

It leverages [optlang](https://github.com/opencobra/optlang) to interface with linear program solvers. This is also the backend that is utilized by [cameo](http://cameo.bio/).

It includes some `matplotlib` plotting functions for generating basic 2D and 3D graphs from the reduced solution spaces (and even 4D if you're willing to use color to represent the fourth dimension).

Some terminology used throught this module:
* **Halfspace or Facet:** A halfspace is a hyper-plane constraint which divides the total solution space in half. In 3D, a halfspace is a bounding plane. In 2D a halfspace is a bounding line. Each halfspace in the final solution is a bounding hyper-plane of the solution space (or a facet).
* **Node:** A node is a unique combination of halfspaces which intersect. These include both facets and vertices depending upon how many halfspaces intersect for a given node. In 3D, a plane node has one halfspace, a line node has two halfspaces, and a vertice node has three halfspaces.
* **Lattice Graph:** This shows the complete set of nodes and the halfspaces they contain

Here's a picture showing a pyramid and the corresponding lattice graph with nodes and halfspaces defined. It also shows how FEA solves the solution; however, read the [full manuscript](docs/_static/FEA_Manuscript.pdf) for more details on how this works.

![Pyramid solution space with lattice graph, facets, and nodes](docs/_static/Pyramid_face_lattice.svg "Pyramid Lattice Graph")