# cut-cell-fem-3D
FEM code based on deal.II (www.dealii.org) to solve 3D problems using unfitted meshes under the cut-cell approach.
The approach used to solve FEM equations in an unfitted mesh is based on the work of Mirtich (1996), where
he outlined a technique to convert integrals on complex polyhedra into line integrals based on Green's theorem.
The classical Poisson problem is solved in a sphere. Nitsche's method was used to impose boundary conditions weakly,
following the mathematical formulation in Burman and Hansbo (2012).

A brief description of the main classes:

NewCell_3D: Contains all the cut cell information, as well as functions to organize and reorder vertices, compute 
geometric center and compute the boundary cut face.

NewFace_3D: Contains all the information about cell lines and faces, both from new cut faces which are generated from the cut
mesh and original faces which maintain deal.ii's native info. It has several methods used to generate the cut faces
and the projected faces, which are used to evaluate relevant integrals.

polymul.h: Modified code from https://code.google.com/p/polymul/. Contains polynomial structure based on dimension
and degree.

Polynomials3D.h: Added several functions to manipulate polynomials, needed to integrate terms on complex polyhedra.

Write_VTK: Creates .vtk files based on lines instead of cells.

References:
Erik Burman and Peter Hansbo. Fictitious domain finite element
methods using cut elements: II. A stabilized Nitsche method. Applied
Numerical Mathematics, 62(4):328 – 341, 2012.

Brian Mirtich. Fast and accurate computation of polyhedral mass
properties. J. Graph. Tools, 1(2):31–50, February 1996.


