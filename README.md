# fem-cut-cell-3D
FEM code based on deal.II (www.dealii.org) to solve 3D problems using unfitted meshes under the cut-cell approach.
The approach used to solve FEM equations in an unfitted mesh is based on the work of Mirtich (1996), where
he outlined a technique to convert integrals on complex polyhedra into line integrals based on Green's theorem.
The classical Poisson problem is solved in a sphere. Nitsche's method was used to impose boundary conditions weakly, following the mathematical formulation in Burman and Hansbo (2012).

The actual code will be released soon after publication. For some results on the 3D Poisson Problem test case, see the folder results.

A brief description of the main classes:

• NewMesh_3D{.cpp,.h}
Responsible for organiing and creating an object of the type dealii::Triangulation<3> containing only elements entirely or partially contained in the main domain.

• NewCell_3D{.cpp,.h}
Contains all the information and all the relevant function for the construction of an intersected element. Holds info about faces through a STL vector of NewFace_3D objects (std::vector<NewFace_3D>) for easy manipulation and data extraction. 

• NewFace_3D{.cpp,.h}

Contains data and functions needed to implement faces of a new cut-cell. Implements several methods used to generate the cut faces and the projected faces, which are used to evaluate relevant integrals. Info about each line are hold in a struct NewLine, which are collected in a STL vector (std::vector<NewLine>).

• Write_VTK {.cpp,.h}

Responsible for generating .vtk files. Generates an structure based in 2D lines, because the newly generated cut-cells may be irregular polyhedra with diverse number of sides.

• polymul.h

Code used to implement and multiply multivariate polynomials, originally found in (Ekstrom, 2009). The code was changed to add some methods to manipulate polynomials that are needed to integrate complex polyhedra.

• Polynomials3D.h

Code containing the class Polynomials3D, which expands the polynomial structure of polymul.h, as well as the NPPolynomials namespace. NPPolynomials contains functions used to manipulate scalar or vector polinomial functions, such as the divergente and gradient operations, computing of vectorial fields, plane projections, matrix-vector multiplication, dot product, etc.

• CutCell_Integration_3D {.cpp,.h}
Implements main functions for numerical integration of terms arising from the finite element formulation using Nitsche's method to impose boundary conditions weakly.



References:
Finite Element formulation:
Erik Burman and Peter Hansbo. Fictitious domain finite element methods using cut elements: II. A stabilized Nitsche method. Applied Numerical Mathematics, 62(4):328 – 341, 2012.

Integration over complex polyhedra:
Brian Mirtich. Fast and accurate computation of polyhedral mass properties. J. Graph. Tools, 1(2):31–50, February 1996.

polymul.h:
EKSTRÖM, U. Polymul library. 2009. Available in: <https://code.google.com/p/polymul/>.



