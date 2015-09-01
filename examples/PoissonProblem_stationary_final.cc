/*
 * Updated with CutCell_Integration.cpp. Reduced total time from 99.8 s to 39.7 s.
 *
 */


#include <ctime>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <math.h> // log()
#include <utility> // pair
#include <vector> // std::veector

#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/base/derivative_form.h>
// If one wants to use ILU preconditioner
#include <deal.II/lac/sparse_ilu.h>

//Added 14/10
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>


#include "../../include/NewCell_3D.h"
#include "../../include/NewMesh_3D.h"
#include "../../include/NewFace_3D.h"
#include "../../include/Write_VTK.h"
#include "../../include/Polynomials3D.h"
#include "../../include/polymul.h"
#include "../../CutCell_Integration_3D.h"

#include <deal.II/dofs/dof_renumbering.h>

namespace cut_cell_method
{
using namespace dealii;
using std::cout;
using std::endl;
enum { X = 0, Y = 1, Z = 2};

//Exact solution of the Poisson problem with zero value B.C.
template <int dim>
class ExactSolution : public Function<dim>
{
public:
	ExactSolution () : Function<dim>() {} // In FAQ:  ExactSolution () : Function<dim>(dim+1) {}

	virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;

	virtual Tensor<1,dim> gradient (const Point<dim>   &p,
			const unsigned int  component = 0) const;
};

template <int dim>
double ExactSolution<dim>:: value (const Point<dim>   &p,
		const unsigned int) const	{
	double return_value = 0;

	return_value = ( 1.0-p.square() )/6.0 ;

	return return_value; // return value = [ ] (size = 1)
}

// Gradient of analytical solution
template <int dim>
Tensor<1,dim > ExactSolution<dim>::gradient (const Point<dim> &p, const unsigned int) const
{
	Tensor<1,dim> return_value;
	return_value = (-p/3);
	return return_value; // return value = [ , ] (size = ?,dim)
}
// Function to "weight" the non-solution elements (elements outside the circle)
template <int dim>
class WeightSolution : public Function<dim>
{
public:
	WeightSolution () : Function<dim>() {} // In FAQ:  ExactSolution () : Function<dim>(dim+1) {}

	virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;
};


template <int dim>
double WeightSolution<dim>:: value (const Point<dim>   &p,
		const unsigned int) const {

	double return_value = 0;
	//	if (p.square() <= Radius_boundary)
	if (p.square() <= 1)
		return_value = 1;
	else
		return_value = 0;

	return return_value; // return value = [ ] (size = 1)
}

//Implementation of a Signed Distance Function of a sphere with radius r.
template<int dim>
class SignedDistanceCircle: public Function<dim> {
	const double radius = 1.0;
	const Point<dim, double> center;
public:
	SignedDistanceCircle() :
		Function<dim>() {
	}
	virtual ~SignedDistanceCircle(){}
	virtual double value(const Point<dim> &p,
			const unsigned int component = 0) const;
};

template<int dim>
double SignedDistanceCircle<dim>::value(const Point<dim> &p,
		const unsigned int component) const {
	Assert(component == 0, ExcInternalError());
	//Center and radius of circle.
	double distance;
	// CIRCLE FUNCTION
	distance = p.distance(center) - radius;

	// ELLIPSE FUNCTION
//	double a_ellipse = 1.5;
//	double b_ellipse = 0.75;
//	distance = (p[0]*p[0]/(a_ellipse*a_ellipse) + p[1]*p[1]/(b_ellipse*b_ellipse)) - 1;

	return distance;
	/*Given any point X0, return (distance - radius) using the following:
	SignedDistanceCircle<dim> example_distance;
				example_distance.value(X0,0)*/
}

// Main class for solution of the Poisson Problem.
template <int dim>
class PoissonProblem
{
public:
	PoissonProblem ();
	void run ();

private:
	enum
	{
		not_surface_domain_id,
		surface_domain_id
	};

	static bool	cell_is_in_bulk_domain (const typename /*hp::*/DoFHandler<dim>::cell_iterator &cell);
	static bool	cell_is_in_surface_domain (const typename /*hp::*/DoFHandler<dim>::cell_iterator &cell);
	void set_active_fe_indices ();

	void make_grid ();
	void initialize_levelset();
	void output_results_levelset() const;
	void get_new_triangulation ();
	void setup_system ();
	void initialize_levelset_new();
	void output_results_levelset_new() const;
	void create_new_cut_cell_mesh();
	void assemble_system();
	void solve ();
	void output_results () const;
	void process_solution();
	void output_results_interpolated();


	int 				 cycle;
	int 				 n_cycles;
	int total_refinement;
	std::string			save_to_folder;
	Triangulation<dim>     triangulation;
	Triangulation<dim>  	 triangulation_new;

	FE_Q<dim>              fe; // <dim> indicates the Dimension
	DoFHandler<dim>        dof_handler;
	DoFHandler<dim>        dof_handler_new;

	MappingQ1 <dim>	 	 mapping;
	SparsityPattern      sparsity_pattern;

	SparseMatrix<double> system_matrix;
	Vector<double>       solution;
	Vector<double>       system_rhs;

	FullMatrix<double>   ALLERRORS;
	Vector<double>       exact_solution;
	Vector<double>       difference_solution;

	Vector<double>       levelset;
	Vector<double>       levelset_project;

	SparsityPattern      sparsity_pattern_new;
	SparseMatrix<double> system_matrix_new;
	SparseMatrix<double> mass_matrix;
	FullMatrix<double>   FM_system_matrix;

	FullMatrix<double> FM_system_matrix_new_with_stab;

	FullMatrix<double>   stabilization_matrix;

	Vector<double>       solution_new;
	Vector<double>       system_rhs_new;

	SparseMatrix<double> system_matrix_real;

	ConstraintMatrix 	 constraints;

	std::vector <bool>   isboundary_face;
	NewMesh_3D 			 Obj_NewMesh;
	std::vector<NewCell_3D> CutTriangulation;
};

template <int dim>
PoissonProblem<dim>::PoissonProblem () // (const FiniteElement<dim> &fe)
:
fe (1), // (1) indicates the polynomial degree
dof_handler ()
, dof_handler_new ()

{}
template <int dim>
bool PoissonProblem<dim>:: cell_is_in_bulk_domain
(const typename /*hp::*/DoFHandler<dim>::cell_iterator &cell)
{
  return (cell->material_id() == not_surface_domain_id);
}


template <int dim>
bool PoissonProblem<dim>::cell_is_in_surface_domain
(const typename /*hp::*/DoFHandler<dim>::cell_iterator &cell)
{
  return (cell->material_id() == surface_domain_id);
}

template <int dim>
void PoissonProblem<dim>::set_active_fe_indices ()
{
	int count_cell = 0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation_new.begin_active();
         cell != triangulation_new.end(); ++cell) {
      if (Obj_NewMesh.vector_new_cell_struct[count_cell].cell_is_boundary) {
        cell->set_material_id (surface_domain_id);
      }
      else
        cell->set_material_id (not_surface_domain_id);
      count_cell++;
    }

    // Set active fe indices only when I upgrade this to the vector-valued problem.
    int inside_cells = 0, surface_cells = 0;
  for (typename /*hp::*/DoFHandler<dim>::active_cell_iterator
       cell = dof_handler_new.begin_active();
       cell != dof_handler_new.end(); ++cell)
    {
      if (cell_is_in_bulk_domain(cell)){
//        cell->set_active_fe_index (1);
        inside_cells++;
      }
      else if (cell_is_in_surface_domain(cell)) {
//        cell->set_active_fe_index (0);
        surface_cells++;
      }
      else
        Assert (false, ExcNotImplemented());
    }
  	std::cout<<" number of cells (inside/surface): "<< inside_cells << " ; "<<surface_cells <<std::endl;
}
template <int dim>
void PoissonProblem<dim>::make_grid ()
{
	std::cout << "Call to make_grid \n";
	if (cycle == 0) // This assures that the Grid is created only once, i.e., in the first cycle;
	{
		GridGenerator::hyper_cube (triangulation, -2, 2);
		triangulation.refine_global(2);// Now, it just makes sense to start in 2
	}
//	int refinement = 0;
//	int refinement = 1;
//	int refinement = 2;
	int refinement = 3;
	total_refinement = 2+refinement;
//	cycle = total_refinement;
	triangulation.refine_global(refinement); // Last refinement cycle
	dof_handler.initialize (triangulation,fe);
	constraints.close();
	dof_handler.distribute_dofs (fe);

	const std::string filename_new = save_to_folder + "/triangulation-ref-"
	+ Utilities::int_to_string(total_refinement, 2) +
	".eps";

	std::ofstream out (filename_new.c_str());
	GridOut grid_out;
	grid_out.write_eps (triangulation, out);
	std::cout << "1st std Triangulation created \n";
	std::cout << "# cells: " << triangulation.n_active_cells() << std::endl;
}

template <int dim>
void PoissonProblem<dim>::initialize_levelset() {
	// levelset function node-wise, ie, the points of the vector are exactly what
	// the function SignedDistanceCircle returns.
	levelset .reinit(dof_handler.n_dofs());
	levelset_project.reinit(dof_handler.n_dofs());

	VectorTools::interpolate(dof_handler, SignedDistanceCircle<dim>(),
			levelset
	);
}

template <int dim>
void PoissonProblem<dim>::output_results_levelset() const {
	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(levelset, "levelset");
	data_out.build_patches();

	const std::string filename = save_to_folder + "/levelset-ref-"
	+ Utilities::int_to_string(total_refinement, 1) +
	".vtk";

	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}
template <int dim>
void PoissonProblem<dim>::get_new_triangulation ()
{
	cout << "Call to get_new_triangulation \n";
	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);

	FEValues<dim> fe_values (fe, quadrature_formula,
			update_values | update_gradients | update_JxW_values
			| update_quadrature_points | update_jacobians |
			update_support_jacobians | update_inverse_jacobians);

	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
			update_values |
			update_gradients |
			update_quadrature_points  |
			update_normal_vectors | update_JxW_values
			/*| update_hessians*/ );


	const unsigned int   dofs_per_cell = fe.dofs_per_cell;
	int count_cell = 0;

	std::vector<types::global_dof_index> cell_global_dof_indices (dofs_per_cell);

	typename
	DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();


	// Creates a vector support_points [n.dofs() x 2] which has the x [0] and y [1]
	// coordinates for each (global) DOF (line)
	std::vector<Point<dim> > support_points(dof_handler.n_dofs());
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);



	bool isinside = false;
	bool isboundary = false;
	Obj_NewMesh.reinit();
	Obj_NewMesh.set_variables(support_points,dofs_per_cell);
	isboundary_face.clear();
	Point<dim> VOID_POINT (-10000,-10000,-10000);
	std::vector< Point<dim> > all_points_vector_firstmesh;
	std::vector<types::global_dof_index> face_dof_indices (GeometryInfo<dim>::vertices_per_face);
	for (; cell!=endc; ++cell)
	{

		fe_values.reinit (cell);
		cell->get_dof_indices(cell_global_dof_indices);

		isinside = false;
		isboundary = false;

		int inside = 0;
		int coincident_node = 0;

		// Identify inside cells (case A)

		for (unsigned int i=0;i<dofs_per_cell;++i) {
			unsigned int j = cell_global_dof_indices[i];

			if (levelset[j]<=0) {
				inside++;
				if (levelset[j] == 0) {
					coincident_node++;
				}
			}
		}
		// Case A: cell is entirely inside boundary
		if (inside == 8)
		{
			isinside = true;
		}
		// Case B: cell is entirely outside boundary, even if 1 node is coincident w/ boundary
		if (inside == 0)
		{
//			isoutside = true;
		}
		if (inside == 1 && coincident_node == 1)
		{
//			isoutside = true;
		}
		// Case C: cell has 3 nodes inside boundary, i.e., is a boundary cell.
		if (inside > 1 && inside < 8 )
		{
			isboundary = true;
		}
		if (inside == 1 && coincident_node != 1 )
		{
			isboundary = true;
		}

		if (isboundary) isboundary_face.push_back(true);
		else isboundary_face.push_back(false);

		bool cell_is_boundary;
		if(isinside)
		{
			cell_is_boundary = false;
			Obj_NewMesh.set_new_mesh(cell_global_dof_indices,
					cell_is_boundary);
		}
		if(isboundary)
		{
			cell_is_boundary = true;
			Obj_NewMesh.set_new_mesh(cell_global_dof_indices,
					cell_is_boundary);
		}
		count_cell++;
	} // (end for cell)

	std::cout << "CYCLE = " << cycle << "\n";
	triangulation_new.clear();
	Obj_NewMesh.create_new_triangulation(triangulation_new);

	const std::string filename_new = save_to_folder + "/triangulation_new-ref-"
	+ Utilities::int_to_string(total_refinement, 1) +
	".eps";

	std::ofstream out (filename_new);
	GridOut grid_out;
	grid_out.write_eps (triangulation_new, out);
	std::cout << "New Triangulation created \n";
	std::cout << "# cells: " << triangulation_new.n_active_cells() << std::endl;
}
template <int dim>
void PoissonProblem<dim>::setup_system ()
{
	set_active_fe_indices();
	dof_handler_new.distribute_dofs (fe);
	DoFRenumbering::Cuthill_McKee (dof_handler_new);
	std::cout << "Number of degrees of freedom: "
			<< dof_handler_new.n_dofs()
			<< std::endl;

	CompressedSparsityPattern c_sparsity_new(dof_handler_new.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler_new, c_sparsity_new);

	sparsity_pattern_new.copy_from(c_sparsity_new);

	system_matrix_new.reinit (sparsity_pattern_new);
	solution_new.reinit (dof_handler_new.n_dofs());
	system_rhs_new.reinit (dof_handler_new.n_dofs());
	FM_system_matrix.reinit (dof_handler_new.n_dofs(),dof_handler_new.n_dofs());

	stabilization_matrix.reinit (dof_handler_new.n_dofs(),dof_handler_new.n_dofs());

	exact_solution.reinit (dof_handler_new.n_dofs());
	difference_solution.reinit (dof_handler_new.n_dofs());

	system_matrix_real.reinit (sparsity_pattern_new);

}

template <int dim>
void PoissonProblem<dim>::initialize_levelset_new() {
	// Need to reinitialize the levelset function with different size, corresponding
	// to the new triangulation.
	levelset.reinit(dof_handler_new.n_dofs());
//	levelset_project.reinit(dof_handler.n_dofs());
	VectorTools::interpolate(dof_handler_new, SignedDistanceCircle<dim>(),
			levelset
			);
}

template <int dim>
void PoissonProblem<dim>::output_results_levelset_new () const {
	DataOut<dim> data_out;
	data_out.attach_dof_handler(dof_handler_new);
	data_out.add_data_vector(levelset, "levelset_new");
	data_out.build_patches();

	const std::string filename = save_to_folder + "/levelset_new-ref-"
			+ Utilities::int_to_string(total_refinement, 1) +
			".vtk";

	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}

template <int dim>
void PoissonProblem<dim>::create_new_cut_cell_mesh()
{
	const clock_t begin_time = clock();
	std::cout << "Call to create_new_cut_cell_mesh \n";
	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);
	QGauss<dim-2> line_quadrature_formula(2);

	FEValues<dim> fe_values (fe, quadrature_formula,
			update_values | update_gradients | update_JxW_values
			| update_quadrature_points /*| update_jacobians | update_support_jacobians
			| update_inverse_jacobians*/);

	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
			update_values |
			update_gradients |
			update_quadrature_points  |
			update_normal_vectors/* | update_JxW_values
			| update_hessians*/ );

	std::vector<types::global_dof_index> cell_global_dof_indices (GeometryInfo<dim>::vertices_per_cell);
	std::vector<types::global_dof_index> face_dof_indices (GeometryInfo<dim>::vertices_per_face);
	std::vector<types::global_dof_index> line_dof_indices (2);
	std::vector<Point<dim> > support_points(dof_handler_new.n_dofs());
	DoFTools::map_dofs_to_support_points(mapping, dof_handler_new, support_points);

	// Create a vector of Objects NewCell, used to retrieve info about cut cells, such as points of
		// intersection, normal vector, length, etc.
	CutTriangulation.reserve(triangulation_new.n_cells());

	std::vector<double> cell_quadrature_weights ;
	cell_quadrature_weights = quadrature_formula.get_weights();

	typename
	/*hp::*/DoFHandler<dim>::active_cell_iterator
	cell = dof_handler_new.begin_active(),
	endc = dof_handler_new.end();
	for (; cell!=endc; ++cell)
	{
		fe_values.reinit(cell);

		cell->get_dof_indices(cell_global_dof_indices);

		if (cell_is_in_bulk_domain(cell))
		{
			NewCell_3D CutCell(/*is_surface_cell*/false);
			CutCell.SetIndex(cell);
			CutTriangulation.push_back(CutCell);
		}

		if (cell_is_in_surface_domain(cell))
		{
			NewCell_3D CutCell(/*is_surface_cell*/true);
			CutCell.SetIndex(cell);

			// Loop over faces, to identify faces
			NewFace_3D CutFaceBoundary(mapping,cell);
			for (unsigned int face_it=0; face_it<GeometryInfo<dim>::faces_per_cell; ++face_it)
			{
				CutCell.InputNewRealFaceInfo(face_it);
				bool face_is_intersected = false;
				bool face_is_inside = false;
				fe_face_values.reinit (cell, face_it);
				cell->face(face_it)->get_dof_indices(face_dof_indices);

				// Loop over lines of the face, if at least one line is crossed by the boundary, this face is intersected.
				for (unsigned int line_it = 0; line_it<GeometryInfo<dim>::lines_per_face; ++line_it)
				{
					cell->face(face_it)->line(line_it)->get_dof_indices(line_dof_indices);
					unsigned int k0 = line_dof_indices[0] ; // Global DOF
					unsigned int k1 = line_dof_indices[1] ;
					if (
							((levelset[k0]*levelset[k1]) < 0.0) ||
						     (levelset[k0] == 0.0 && levelset[k1] < 0.0) ||
													(levelset[k1] == 0.0 && levelset[k0] < 0.0) )
					{
						face_is_intersected = true;
					}
				}
//				 Loop over lines of a no intersected face, to retrieve inside faces (inside faces with coincident nodes should be counted as intersected faces)
				if (!face_is_intersected)
					for (unsigned int line_it = 0; line_it<GeometryInfo<dim>::lines_per_face; ++line_it)
					{
						cell->face(face_it)->line(line_it)->get_dof_indices(line_dof_indices);
						unsigned int k0 = line_dof_indices[0] ; // Global DOF
						unsigned int k1 = line_dof_indices[1] ;
						if (	(levelset[k0] < 0.0 )
								&& (levelset[k1]) < 0.0 )
							face_is_inside = true;
					}

				if (face_is_intersected)
				{

					std::vector<Point<dim> > new_boundary_CutLine(2);
					NewFace_3D CutFace(levelset,mapping,cell);
					CutCell.DefStabilizationFace (face_it, true);
					CutFace.SetStabilizationFace(true); // Set this as a face that will be later integrated in the stabilization term.
					int integrate_line_check = 0;  // keep track of the number of the integrated
					int intersected_face = 0; // keep track of the number of intersected faces;
					for (unsigned int line_it = 0; line_it<GeometryInfo<dim>::lines_per_face; ++line_it)
					{
						cell->face(face_it)->line(line_it)->get_dof_indices(line_dof_indices);

//						faces of the cell, i.e., new cut faces and other boundary faces of the cell

						int inside_faces = 0; // keep track of the number of inside faces;
						unsigned int k0 = line_dof_indices[0] ; // Global DOF
						unsigned int k1 = line_dof_indices[1] ;

						if ((levelset[k0]*levelset[k1]) < 0.0)
						{
							Point<dim> X0;
							Point<dim> X1;
							Point<dim> face_normal_vector;

							X1 = CutFace.GetIntersection(support_points[k0],support_points[k1], k0,k1);

							if (levelset[k1] > levelset[k0])
								X0 = support_points[k0];
							else X0 = support_points[k1];
							// Setting the points of the new face.
							new_boundary_CutLine[intersected_face] = X1;
							intersected_face++;

							// This is to prevent the integration when the boundary
							// reaches a DOF such that a new face is created with
							// the same two DOF's.
							if (X0 != X1)
							{
								// This function also Inputs new Line
								CutFace.SetCoordinates(integrate_line_check,X0,X1,false,true);
								CutFace.SetUnitCoordinates(integrate_line_check,
								mapping.transform_real_to_unit_cell(cell,	X0),
								mapping.transform_real_to_unit_cell(cell,	X1));
								integrate_line_check++;
							}
						} // End if line is intersected by the boundary.

						// The next 2 statements are only to include the coincident node as
						// a node for the new face. (It is in fact the intersection)
						if (fabs(levelset[k0]) < pow(10,-10) && levelset[k1] < 0.0)
						{
//							if (levelset[k0] == 0.0 && levelset[k1] < 0.0) {
							new_boundary_CutLine[intersected_face] = support_points[k0];
							intersected_face++;
						}
						else if (fabs(levelset[k1]) < pow(10,-10)  && levelset[k0] < 0.0) {
//							else if (levelset[k1] ==0  && levelset[k0] < 0.0) {
//							std::cout << "--Line is inside the domain with one node coincident (k1)\n";
							new_boundary_CutLine[intersected_face] = support_points[k1];
							intersected_face++;
						}

//						Find INSIDE FACES, not intersected by the boundary (but inside bulk domain).
//						Exclude faces with both DOFs outside the boundary.
						if (levelset[k0]<=0.0 && levelset[k1]<=0.0)
						{
							Point<dim> X0;
							Point<dim> X1;
							Point<dim> face_normal_vector;
							inside_faces++; // Keep track of number of inside_faces.
							X0 = support_points[k0];
							X1 = support_points[k1];
							// This function also Inputs new Line
							CutFace.SetCoordinates(integrate_line_check,X0,X1,false,true);
							CutFace.SetUnitCoordinates(integrate_line_check,
															mapping.transform_real_to_unit_cell(cell,	X0),
															mapping.transform_real_to_unit_cell(cell,	X1));
							integrate_line_check++;
						} // End Find INSIDE Faces.
					} // end loop lines
					// Create NEW CUT FACE (BOUNDARY)
					Point<dim> X0;
					Point<dim> X1;
					Point<dim> face_normal_vector;

					X0 = new_boundary_CutLine[0];
					X1 = new_boundary_CutLine[1];
					CutFaceBoundary.Add1Line();
					// This function also Inputs new Line
					CutFace.SetCoordinates(integrate_line_check,X0,X1,true,true);
					CutFace.SetUnitCoordinates(integrate_line_check,
													mapping.transform_real_to_unit_cell(cell,	X0),
													mapping.transform_real_to_unit_cell(cell,	X1));

					// All lines of the face were set; now need to compute the centroid of the face.
					CutFace.CompFaceCentroid();
					CutFace.SetFaceNormal(fe_face_values.normal_vector(0));
//					CutFace.CompProjectionVars();
					CutCell.InputNewFace(CutFace);
				} // end if face is intersected
//				If face is not intersected, but still is part of the cut cell, it needs to be computed!
//				However, I only want INSIDE FACES.
				else if (face_is_inside)
				{
					NewFace_3D CutFace(levelset,mapping,cell);
					CutFace.SetStabilizationFace(true);
					CutCell.DefStabilizationFace (face_it, true);
					int integrate_line_check = 0;
					for (unsigned int line_it = 0; line_it<GeometryInfo<dim>::lines_per_face; ++line_it)
					{
						cell->face(face_it)->line(line_it)->get_dof_indices(line_dof_indices);
						unsigned int k0 = line_dof_indices[0] ; // Global DOF
						unsigned int k1 = line_dof_indices[1] ;

						Point<dim> X0 = support_points[k0];
						Point<dim> X1 = support_points[k1];
						Point<dim> face_normal_vector;

						CutFace.SetCoordinates(integrate_line_check,X0,X1,/*is_boundary_face*/false,/*add_1_line*/true);
						CutFace.SetUnitCoordinates(integrate_line_check,
														mapping.transform_real_to_unit_cell(cell,	X0),
														mapping.transform_real_to_unit_cell(cell,	X1));
						++integrate_line_check;
					} // end for line
					// All lines of the face were set; now need to compute the centroid of the face.
					CutFace.CompFaceCentroid();
					CutFace.SetFaceNormal(fe_face_values.normal_vector(0));
//					CutFace.CompProjectionVars();
					CutCell.InputNewFace(CutFace);
				} // end if face is inside
			   // end if face is not intersected
			} // end loop faces
//			Found all intersections and set all new cut faces; Now need to create the boundary cut face.
			assert(CutCell.real_faces_vector.size() == GeometryInfo<dim>::faces_per_cell && "Real faces were not input correctly");
			CutCell.CompBoundaryFace(CutFaceBoundary); // Call to SetCoordinates; Call to SetVertices
			for (unsigned int line_it = 0;
					line_it < /*CutCell.Obj_VectorNewFace[CutCell.number_of_faces-1]*/CutFaceBoundary.number_of_lines ; ++line_it)
			{

				CutFaceBoundary.SetUnitCoordinates // This also calls SetUnitVertices
				(		line_it,
						mapping.transform_real_to_unit_cell(cell, CutFaceBoundary.Obj_VectorNewLine[line_it].X0 ),
						mapping.transform_real_to_unit_cell(cell, CutFaceBoundary.Obj_VectorNewLine[line_it].X1 )
						);
			} // end for lines


			CutFaceBoundary.CompFaceCentroid(); // This also calls CompUnitFaceCentroid
			CutFaceBoundary.SetStabilizationFace(false);
			// Only now that I have all the info about the faces I can calculate the normal vectors - the calculation needs the Centroid of the face.

			//			Compute the cell_centroid of the cell, needed to organize the vertices of each line and the normal vector of each face.
						CutCell.CompCellCentroid(); // Calls OrganizeVertices
			CutFaceBoundary.CompCutFaceNormal(CutCell.cell_centroid);
			CutCell.InputNewFace(CutFaceBoundary);
			for (unsigned int face_it = 0; face_it < CutCell.number_of_faces ; ++face_it)
			{
				CutCell.Obj_VectorNewFace[face_it].CompAllLineNormals();
				CutCell.Obj_VectorNewFace[face_it].CompProjectionVars();
			}

			assert(CutCell.Obj_VectorNewFace[CutCell.number_of_faces-1].is_boundary_face == true)	;

//			Set the vertices of each cell so they are in a counter clockwise order. Need the normal vector of the face, therefore need to be
//			called AFTER CompCutFaceNormal
			CutCell.ReorderAllVertices();
			// Calculate the VOLUME of each cell.
			CutCell.CompCellVolume();
			CutTriangulation.push_back(CutCell);
		} // end cell is in surface domain
	} // end for cell
	std::cout << "End loop cells \n";

//	Create a triangulation with stabilization faces only, for debugging purposes.
	if(0)
	{
		std::vector<NewCell_3D> CutTriangulation_only_stab_faces;
		int cell_it;
		for (cell = dof_handler_new.begin_active(), cell_it = 0; cell!=endc; ++cell, ++cell_it)
		{
			if (CutTriangulation[cell_it].is_surface_cell)
			{
				NewCell_3D CutCell(/*is_surface_cell*/true);
				CutCell.SetIndex(cell);

				for (int face_it = 0; face_it < GeometryInfo<dim>::faces_per_cell; ++face_it)
				{
					if (CutTriangulation[cell_it].real_faces_vector[face_it].is_stabilization_face == true)
					{
						NewFace_3D CutFace(levelset,mapping,cell);
						int integrate_line_check = 0;
						for (unsigned int line_it = 0; line_it<GeometryInfo<dim>::lines_per_face; ++line_it)
						{
							cell->face(face_it)->line(line_it)->get_dof_indices(line_dof_indices);

							unsigned int k0 = line_dof_indices[0] ; // Global DOF
							unsigned int k1 = line_dof_indices[1] ;

							Point<dim> X0;
							Point<dim> X1;
							//						inside_faces++; // Keep track of number of inside_faces.
							X0 = support_points[k0];
							X1 = support_points[k1];
							CutFace.SetCoordinates(integrate_line_check,X0,X1,false,true);
							++integrate_line_check;
						} // end loop lines.
						CutCell.InputNewFace(CutFace);
					} // end if face is stab face
				} // end loop faces
				CutTriangulation_only_stab_faces.push_back(CutCell);
			} // end if is surface cell
		} // end loop cell
	}



//	Output CutTriangulation as VTK file.
	{
		Write_VTK Obj_Write_VTK(CutTriangulation);
	Obj_Write_VTK.OrganizeVertices();

	std::string filename_new = save_to_folder + "/Cut-cell_mesh-ref-";
	filename_new += ('0' + total_refinement);
	filename_new += ".vtk";
	std::ofstream output(filename_new.c_str());
}

	if(0) // Create a triangulation with only stabilization faces to check if the faces were assigned correctly. Everything seems in order.
	 {
		 Write_VTK Obj_Write_VTK_OneCell(CutTriangulation_only_stab_faces);
		 Obj_Write_VTK_OneCell.OrganizeVertices();
		 std::string filename_new = save_to_folder + "/Cut-cell_mesh_only_stab_faces-ref-";
		 filename_new += ('0' + total_refinement);
		 filename_new += ".vtk";
		 std::ofstream output_OneCell(filename_new.c_str());

		 Obj_Write_VTK_OneCell.WriteFile(output_OneCell);
		 Obj_Write_VTK_OneCell.AddVectors(output_OneCell);
		 cout << "CutTriangulation_only_stab_faces created \n";
	 }

 std::cout << "Call to create_new_cut_cell_mesh successfull.\n";
 cout << "Time elapsed: \n";
 std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
 cout << endl;
}

template <int dim>
void PoissonProblem<dim>::assemble_system()
{
	std::cout << "Call to assemble_system \n";
	const clock_t begin_time = clock();
	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-1> face_quadrature_formula(2);
	QGauss<dim-2> line_quadrature_formula(2);

	FEValues<dim> fe_values (fe, quadrature_formula,
			update_values | update_gradients | update_JxW_values
			| update_quadrature_points | update_jacobians | update_support_jacobians
			| update_inverse_jacobians);

	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
			update_values |
			update_gradients |
			update_quadrature_points  |
			update_normal_vectors | update_JxW_values	);

    FEFaceValues<dim>   fe_face_values_neighbor_cell (fe,
            											face_quadrature_formula,
                                                        update_JxW_values |
                                                        update_normal_vectors |
                                                        update_gradients);

	// Creates a vector support_points [n.dofs() x 2] which has the x [0] and y [1]
		// coordinates for each (global) DOF (line)
//		std::vector<Point<dim> > support_points(dof_handler.n_dofs());
//		DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);
		const unsigned int   dofs_per_cell = fe.dofs_per_cell;
		const unsigned int   n_q_points    = quadrature_formula.size();
		const unsigned int   n_face_q_points    = face_quadrature_formula.size();
		const unsigned int   line_n_q_points = line_quadrature_formula.size();
		Vector<double>       cell_rhs (dofs_per_cell);
		FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
		FullMatrix<double>   cell_stabilization_matrix (dofs_per_cell+4, dofs_per_cell+4);
		std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
				std::vector<types::global_dof_index> face_dof_indices (GeometryInfo<dim>::vertices_per_face);

				std::vector<double> face_quadrature_weights;
				face_quadrature_weights = face_quadrature_formula.get_weights();

		std::vector<double> line_quadrature_weights;
		line_quadrature_weights = line_quadrature_formula.get_weights();
		std::vector<Point<1> > line_quadrature_points;
		line_quadrature_points = line_quadrature_formula.get_points();
	//Test Poi evaluating poly phi_0.
	FullMatrix<double>  jacobian_inverse(dim,dim);
//	Nitsche Method parameter.
	const double gamma_D = 5.0;
//	Stabilization parameter.
	const double gamma_1 = 0.1;
//	RHS Value.
	const double f_B = 1.0;
	double cell_diameter;

	std::vector<double> cell_quadrature_weights ;
	cell_quadrature_weights = quadrature_formula.get_weights();
	std::vector<types::global_dof_index> line_dof_indices (2);
	double jacobian_determinant ; /*= fe_values.JxW(0)/(cell_quadrature_weights[0])*/
	std::vector<types::global_dof_index> cell_global_dof_indices (GeometryInfo<dim>::vertices_per_cell);

	assert (line_n_q_points == 2);

		std::vector<Point<dim> > support_points(dof_handler_new.n_dofs());
		DoFTools::map_dofs_to_support_points(mapping, dof_handler_new, support_points);
		std::vector<types::global_dof_index> neighbor_cell_global_dof_indices /*(dofs_per_cell)*/;
		std::vector<types::global_dof_index> neighbor_face_global_dof_indices /*(dofs_per_cell/2)*/;
		std::vector<int> j_integrate_face;

		CutCell_Integration_3D Obj_CutCell_Integration_3D (CutTriangulation, dofs_per_cell);

		typename DoFHandler<dim>::active_cell_iterator
							cell = dof_handler_new.begin_active(),	endc = dof_handler_new.end();
		int cell_it = 0;
		for (cell = dof_handler_new.begin_active(),  cell_it = 0; cell!=endc; ++cell,++cell_it)
		{
			cout << "Cell: " << cell << " cell_it: " << cell_it << endl;
//			CutTriangulation[cell_it].OutputVertices();
			fe_values.reinit (cell);
			jacobian_determinant = fe_values.JxW(0)/(cell_quadrature_weights[0]);
			cell_diameter = cell->diameter();

			jacobian_inverse = 0;
			for (unsigned int i=0; i<dim; ++i)
				for (unsigned int j=0; j<dim; ++j)
					jacobian_inverse(i,j) = fe_values.inverse_jacobian(0)[i][j];

			cell_matrix = 0.0;
			cell_rhs = 0.0;

			cout <<  "NEW CELL ---------------------------------------------------  " << cell << endl;
			if (CutTriangulation[cell_it].is_surface_cell)
			{
				assert (cell_is_in_surface_domain(cell));
				Obj_CutCell_Integration_3D.initialize_cell_variables(cell_it, jacobian_determinant, jacobian_inverse);

				for (int face_it = 0; face_it < CutTriangulation[cell_it].number_of_faces; ++face_it)
				{
					Obj_CutCell_Integration_3D.initialize_face_variables (face_it);

					for (int line_it = 0; line_it <  CutTriangulation[cell_it].Obj_VectorNewFace[face_it].number_of_lines; ++line_it)
					{
//
						Obj_CutCell_Integration_3D.initialize_line_variables (line_it);
						for (unsigned int dof_i = 0; dof_i<dofs_per_cell; ++dof_i)
						{
							for (unsigned int dof_j = 0; dof_j<dofs_per_cell; ++dof_j)
							{
//								Integrate bilinear term a(phi_i,phi_j)
								cell_matrix(dof_i,dof_j) += Obj_CutCell_Integration_3D.compute_bilinear_polynomial(dof_i,dof_j);
//								Integrate boundary terms of Nitsche's formulation (terms C, D, D2)
								if (CutTriangulation[cell_it].Obj_VectorNewFace[face_it].is_boundary_face)
									cell_matrix(dof_i,dof_j) +=Obj_CutCell_Integration_3D.compute_terms_c_d(dof_i,dof_j,gamma_D,
											cell_diameter);
							}
//							Integrate term b of RHS = (f,phi) dOmega.
							if (f_B != 0.0)
							{
								cell_rhs(dof_i) += Obj_CutCell_Integration_3D.compute_rhs(dof_i,f_B);
							}

						}
					} // End loop lines
				} // End loop faces
				// Integrate stabilization terms. j(u,v)
				// Need to loop over REAL faces, not cut faces.
				for (unsigned int face_it = 0; face_it < GeometryInfo<dim>::faces_per_cell ; ++face_it)
				{
					if (CutTriangulation[cell_it].real_faces_vector[face_it].is_stabilization_face == true)
					{
					fe_face_values.reinit(cell,face_it);
					cell->face(face_it)->get_dof_indices(face_dof_indices);
//					Loop over neighbor cells.
					for (unsigned int neighbor_cell_it =0; neighbor_cell_it < GeometryInfo<dim>::faces_per_cell; ++neighbor_cell_it)
					{
//						Loop over 6 faces, if neighbor idx is -1, it doesn't exist (main cell is at the boundary).
						if (cell->neighbor_index(neighbor_cell_it) != -1)
						{
//							Loop over neighbor faces.
							for (unsigned int neighbor_cell_face_it = 0; neighbor_cell_face_it < GeometryInfo<dim> ::faces_per_cell;
									++neighbor_cell_face_it)
							{
//								If neighbor cell's face is common with main cell's face.
								if (cell->neighbor(neighbor_cell_it) // neighbor cell K'
										->face(neighbor_cell_face_it) // faces of K' cell
										== cell->face(face_it) )
									// This means that the cell K' neighbor_cell_it has a
									// face in common with the present cell K.
								{
									fe_face_values_neighbor_cell.reinit(cell->neighbor(neighbor_cell_it),
											neighbor_cell_face_it);

									int neighbor_dofs_per_cell
									= cell->neighbor(neighbor_cell_it)->get_fe().dofs_per_cell;

									int neighbor_dofs_per_face
									= cell->neighbor(neighbor_cell_it)->get_fe().dofs_per_face;

									neighbor_cell_global_dof_indices.resize(neighbor_dofs_per_cell);
									neighbor_face_global_dof_indices.resize(neighbor_dofs_per_face);
									cell->get_dof_indices(cell_global_dof_indices);
									// Get DOF indices of the neighbor CELL
									cell->neighbor(neighbor_cell_it) // neighbor cell K'
												->get_dof_indices(neighbor_cell_global_dof_indices);

									// Get DOF indices of the neighbor FACE (is common to the K face)
									cell->neighbor(neighbor_cell_it)
												->face(neighbor_cell_face_it)->get_dof_indices(neighbor_face_global_dof_indices);
									int face_global_index = cell->face_index(face_it);
//									 Find if this face has already been submitted to the stabilization term integration.
//									j_integrate_face is a vector storing global indices of each face, so that I can keep track of the faces that
//									have already been integrated.
									// If this is true, face has not been found.
									if (std::find(j_integrate_face.begin(),
											j_integrate_face.end(), face_global_index)
									== j_integrate_face.end())
									{
//											cout << "Face was not J-integrated yet \n";
										j_integrate_face.push_back(face_global_index);

										// Create a vector EXTENDED_global_dof_indices with the same global dof's as cell K...
										std::vector<int> EXTENDED_global_dof_indices;
										for (unsigned int i=0;i<dofs_per_cell;++i)
										{
											EXTENDED_global_dof_indices.push_back(cell_global_dof_indices[i]);
//												cout << "Push_back global dof index: " << cell_global_dof_indices[i] << endl;
//												cout << "EXTENDED_global_dof_indices.size(): " << EXTENDED_global_dof_indices.size() << endl;
										}

//										 Check which DOF's of cell K' don't belong to cell K and add them to EXTENDED_global_dof_indices.
//										EXTENDED_global_dof_indices now will have all the DOF's of K and the DOF's of its neighbor K'.
										for (unsigned int i=0;i<dofs_per_cell;++i)
										{
											if (std::find(EXTENDED_global_dof_indices.begin(),
													EXTENDED_global_dof_indices.end(),
													neighbor_cell_global_dof_indices[i])
											== EXTENDED_global_dof_indices.end() )
											{
												EXTENDED_global_dof_indices.push_back( neighbor_cell_global_dof_indices[i] );
//													cout << "Push_back (neighbor) global dof index: " << neighbor_cell_global_dof_indices[i] << endl;
//													cout << "EXTENDED_global_dof_indices.size(): " << EXTENDED_global_dof_indices.size() << endl;
											}
										}

//											cout << "EXTENDED_global_dof_indices.size(): " << EXTENDED_global_dof_indices.size() << endl;
										assert(EXTENDED_global_dof_indices.size() == (dofs_per_cell+dofs_per_cell/2));

//										Assert that the same first global indices of K are the same as K' 's.
										assert (face_dof_indices[0] == neighbor_face_global_dof_indices[0]);
										assert (face_dof_indices[1] == neighbor_face_global_dof_indices[1]);
										assert (face_dof_indices[2] == neighbor_face_global_dof_indices[2]);
										assert (face_dof_indices[3] == neighbor_face_global_dof_indices[3]);


										// Select the local DOF indices of the face in the neighbor cell K'
										// corresponding to the equivalent global DOF indices in the
										// face in cell cell K.
										std::vector<int> EXTENDED_local_dof_K;
										std::vector<int> EXTENDED_local_dof_K_neighbor;

										for (unsigned int i=0;i<dofs_per_cell;++i)
											EXTENDED_local_dof_K.push_back(i);

										EXTENDED_local_dof_K.push_back(-1); // 9
										EXTENDED_local_dof_K.push_back(-1);	// 10
										EXTENDED_local_dof_K.push_back(-1);	// 11
										EXTENDED_local_dof_K.push_back(-1);	// 12

										for (unsigned int i=0;i<EXTENDED_global_dof_indices.size() ;++i)
										{
											// Is global index (i) of neighbor cell K' a
											// global dof of cell K? If NOT:
											if (std::find(neighbor_cell_global_dof_indices.begin(),
													neighbor_cell_global_dof_indices.end(),
													EXTENDED_global_dof_indices[i])
											== neighbor_cell_global_dof_indices.end())
											{
												EXTENDED_local_dof_K_neighbor.push_back(-1);
											}
											else
											{
												for (unsigned int j=0;j<dofs_per_cell;++j)
												{
													if (neighbor_cell_global_dof_indices[j] ==
															EXTENDED_global_dof_indices[i])
														EXTENDED_local_dof_K_neighbor.push_back(j);
												}
											}
										}
										// Important, never forget to reset local element matrix!
										cell_stabilization_matrix = 0.0;
										for (unsigned int dof_i = 0; dof_i<dofs_per_cell+4; ++dof_i)
										{
											for (unsigned int dof_j = 0; dof_j<dofs_per_cell+4; ++dof_j)
											{
													double j = 0;
													double dof_i_jump = 0;
													double dof_j_jump = 0;

													for (int q_point = 0;q_point<n_face_q_points;++q_point)
													{
														Point<dim> normal = fe_face_values.normal_vector(q_point);
														// Commom global DOF's for K and K'
														if ( EXTENDED_local_dof_K[dof_i] != -1 && EXTENDED_local_dof_K_neighbor[dof_i] != -1)
														{
															dof_i_jump =  normal *fe_face_values.shape_grad(EXTENDED_local_dof_K[dof_i] ,q_point) -
																	normal * fe_face_values_neighbor_cell.shape_grad(EXTENDED_local_dof_K_neighbor[dof_i] ,q_point);
														}
														// Global DOF's are not in the cell K
														else if ( EXTENDED_local_dof_K[dof_i] == -1 && EXTENDED_local_dof_K_neighbor[dof_i] != -1)
														{
															dof_i_jump = ( 0 -
																	normal * fe_face_values_neighbor_cell.shape_grad(EXTENDED_local_dof_K_neighbor[dof_i] ,q_point)) ;
														}
														// Global DOF's are not in the cell K'
														else if ( EXTENDED_local_dof_K_neighbor[dof_i] == -1 && EXTENDED_local_dof_K[dof_i] != -1)
														{
															dof_i_jump = ( normal *fe_face_values.shape_grad(EXTENDED_local_dof_K[dof_i] ,q_point) - 0) ;
														}
														// No common global DOF - should not happen
														else assert(0 && "No common global DOF, should not happen");

														// Common global DOF's for K and K'
														if ( EXTENDED_local_dof_K[dof_j] != -1 && EXTENDED_local_dof_K_neighbor[dof_j] != -1)
														{
															dof_j_jump =  normal *fe_face_values.shape_grad(EXTENDED_local_dof_K[dof_j] ,q_point) -
																	normal * fe_face_values_neighbor_cell.shape_grad(EXTENDED_local_dof_K_neighbor[dof_j] ,q_point);
														}
														// Global DOF's are not in the cell K
														else if ( EXTENDED_local_dof_K[dof_j] == -1 && EXTENDED_local_dof_K_neighbor[dof_j] != -1)
														{
															dof_j_jump = ( 0 -
																	normal * fe_face_values_neighbor_cell.shape_grad(EXTENDED_local_dof_K_neighbor[dof_j] ,q_point) );
														}
														// Global DOF's are not in the cell K'
														else if ( EXTENDED_local_dof_K_neighbor[dof_j] == -1 && EXTENDED_local_dof_K[dof_j] != -1)
														{
															dof_j_jump = ( normal *fe_face_values.shape_grad(EXTENDED_local_dof_K[dof_j] ,q_point) - 0) ;
														}
														else assert(0 &&  "No common global DOF, should not happen");
														j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);

													}  // End quadrature points
													cell_stabilization_matrix(dof_i, dof_j )+= j/**gamma_1*cell_diameter*/;
//													gamma_1*cell_diameter is now in the sum of system_matrix+stab. matrix
											} // end dof_j
										} // end dof_i
										for (unsigned int i=0; i<dofs_per_cell+4; ++i)
											for (unsigned int j=0; j<dofs_per_cell+4; ++j)
												stabilization_matrix.add(EXTENDED_global_dof_indices[i],
														EXTENDED_global_dof_indices[j],
														cell_stabilization_matrix(i,j));
									} // end if face was found (I mean, if face was already j-integrated)
//										else cout << "Face was ALREADy J-integrated \n";
								} // end if face is common
							} // end loop faces of neighbor cell.
						}  // end if neighbor cell is !=-1
					} // end loop neighbor cells
				} // end if face is part of stabilization term
			} // end loop real faces
			} // End if cell is in surface domain
			else if (CutTriangulation[cell_it].is_surface_cell == false)
			{
				assert (cell_is_in_bulk_domain(cell));
				fe_values.reinit (cell);
				for (unsigned int q_index = 0; q_index< n_q_points ; ++q_index)
					for (unsigned int dof_i = 0; dof_i<dofs_per_cell; ++dof_i)
					{
						for (unsigned int dof_j=0; dof_j<dofs_per_cell; ++dof_j)
							cell_matrix(dof_i,dof_j) += fe_values.shape_grad(dof_i,q_index)*fe_values.shape_grad(dof_j,q_index)
							*fe_values.JxW (q_index);

						cell_rhs(dof_i) += f_B*fe_values.shape_value (dof_i, q_index)*fe_values.JxW (q_index);
					}
			}
			else assert(0 && "Cell is neither surface nor bulk; something's wrong");

				cell->get_dof_indices (local_dof_indices);
				for (unsigned int dof_i = 0; dof_i<dofs_per_cell; ++dof_i)
			      for (unsigned int dof_j=0; dof_j<dofs_per_cell; ++dof_j)
			    	  system_matrix_new.add (local_dof_indices[dof_i],
			                           local_dof_indices[dof_j],
			                           cell_matrix(dof_i,dof_j));

				for (unsigned int dof_i = 0; dof_i<dofs_per_cell; ++dof_i)
					system_rhs_new(local_dof_indices[dof_i]) += cell_rhs(dof_i);

				ExactSolution<dim> extract_exact_solution;
				//Creating the vector of exact solution:
				for (unsigned int i=0; i<dofs_per_cell; ++i)
				{
					unsigned int j = local_dof_indices[i];
					exact_solution[j] =  extract_exact_solution.value(support_points[j]);
				}

		} // end for cell (DofHandler)

				ALLERRORS[cycle][3] = cell_diameter;
				cout << "gamma_D: " << gamma_D << " gamma_1: " << gamma_1 << endl;
				cout << "cell_diameter:  " << cell_diameter << endl;

//				Calculate L1 condition number of the stiffness matrix, without stabilization.
				{
					FullMatrix<double> IFM_system_matrix(dof_handler_new.n_dofs(),dof_handler_new.n_dofs());
					FullMatrix<double> FM_system_matrix_new(dof_handler_new.n_dofs(),dof_handler_new.n_dofs());
					FM_system_matrix_new.copy_from(system_matrix_new);
					IFM_system_matrix.copy_from(FM_system_matrix_new);
					IFM_system_matrix.gauss_jordan();
					double condition_number = IFM_system_matrix.l1_norm()*FM_system_matrix_new.l1_norm();
					std::cout << "L1 Condition Number without stabilization: " <<  condition_number << "\n";
					//				ALLERRORS[cycle][7] = condition_number;
				}

//				Create system_matrix with stabilization terms.
				FM_system_matrix_new_with_stab.reinit(dof_handler_new.n_dofs(),dof_handler_new.n_dofs());
				FM_system_matrix_new_with_stab.copy_from(system_matrix_new);
				FM_system_matrix_new_with_stab.add(gamma_1*cell_diameter,stabilization_matrix);

//				Calculate L1 condition number of the stiffness matrix, with stabilization.
				{
					double condition_number_w ;
					FullMatrix<double> IFM_system_matrix_new_with_stab(dof_handler_new.n_dofs(),dof_handler_new.n_dofs());
					IFM_system_matrix_new_with_stab.copy_from(FM_system_matrix_new_with_stab);
					IFM_system_matrix_new_with_stab.gauss_jordan();
					/*double*/ condition_number_w = IFM_system_matrix_new_with_stab.l1_norm()*FM_system_matrix_new_with_stab.l1_norm();
					std::cout << "L1 Condition Number with stabilization: " <<  condition_number_w << "\n";
					ALLERRORS[cycle][7] = condition_number_w;
				}

			// Output & Visualize the stabilization matrix
				if(0)
				{
					std::ofstream STABILIZATION_MATRIX;
					STABILIZATION_MATRIX.open("stabilization_matrix.txt");
					for(int unsigned i = 0; i < stabilization_matrix.size(0); ++i) {
						for(int unsigned j = 0; j < stabilization_matrix.size(1); ++j) {
							STABILIZATION_MATRIX << stabilization_matrix(i,j) << ',';
						}
						STABILIZATION_MATRIX << std::endl;
					}
					STABILIZATION_MATRIX.close();
				}

			// Output & Visualize the SYSTEM matrix
				if(0)
				{
					FullMatrix<double> FM_system_matrix_new(dof_handler_new.n_dofs(),dof_handler_new.n_dofs());
					FM_system_matrix_new.copy_from(system_matrix_new);
					std::ofstream SYSTEM_MATRIX_NEW;
					SYSTEM_MATRIX_NEW.open("FM_system_matrix_new.txt");
					for(int unsigned i = 0; i < FM_system_matrix_new.size(0); ++i) {
						for(int unsigned j = 0; j < FM_system_matrix_new.size(1); ++j) {
							SYSTEM_MATRIX_NEW << FM_system_matrix_new(i,j) << ',';
						}
						SYSTEM_MATRIX_NEW << std::endl;
					}
					SYSTEM_MATRIX_NEW.close();
				}
				// Output & Visualize the Stabilized stiffness matrix
				if(0)
			{
				std::ofstream FM_SYSTEM_MATRIX_NEW_WITH_STAB;
				std::string filename_new = save_to_folder + "/FM_system_matrix_new_with_stab.txt";
				FM_SYSTEM_MATRIX_NEW_WITH_STAB.open(filename_new.c_str());
				for(int unsigned i = 0; i < FM_system_matrix_new_with_stab.size(0); ++i) {
					for(int unsigned j = 0; j < FM_system_matrix_new_with_stab.size(1); ++j) {
						FM_SYSTEM_MATRIX_NEW_WITH_STAB << FM_system_matrix_new_with_stab(i,j) << ',';
					}
					FM_SYSTEM_MATRIX_NEW_WITH_STAB << std::endl;
				}
				FM_SYSTEM_MATRIX_NEW_WITH_STAB.close();
			}
			cout << "Time elapsed: \n";
			 std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
			 cout << endl;
}

template <int dim>
void PoissonProblem<dim>::solve ()
{
	std::cout << "Solving... \n";

		SolverControl           solver_control (10000, 1e-12);
	SolverCG<>              solver (solver_control);

//	solver.solve (system_matrix_new , solution_new, system_rhs_new,
//			PreconditionIdentity());
	solver.solve (FM_system_matrix_new_with_stab , solution_new, system_rhs_new,
				PreconditionIdentity());
	  std::cout << "   " << solver_control.last_step()
	            << " CG iterations needed to obtain convergence."
	            << std::endl;

}

template <int dim>
void PoissonProblem<dim>::output_results () const
{
	DataOut<dim> data_out_new;
	data_out_new.attach_dof_handler (dof_handler_new);
	data_out_new.add_data_vector (solution_new, "solution_new");
	data_out_new.build_patches ();

	std::string filename_new = save_to_folder + "/solution_new-ref-";
	filename_new += ('0' + total_refinement);
	filename_new += ".vtk";
	std::ofstream output_new(filename_new.c_str());
	data_out_new.write_vtk(output_new);

	{
		// 3D Output of EXACT SOLUTION, one for each cycle of global refinement
		DataOut<dim> data_out_exact;
		data_out_exact.attach_dof_handler (dof_handler_new);
		data_out_exact.add_data_vector (exact_solution, "exact_solution");
		data_out_exact.build_patches ();

		std::string filename_exact = save_to_folder + "/exact_solution-ref-";
		filename_exact += ('0' + total_refinement);
		filename_exact += ".vtk";
		std::ofstream output_exact(filename_exact.c_str());
		data_out_exact.write_vtk(output_exact);
	}
}

template <int dim>
void PoissonProblem<dim>::process_solution ()
{

	cout << "Call to process_solution \n";
	WeightSolution<dim> weight_values;

	const unsigned int n_active_cells=triangulation_new.n_active_cells();
	const unsigned int n_dofs=dof_handler_new.n_dofs();

	Vector<double> weight_integrate_difference(n_active_cells); // Vector initiliazed, all elements equal to zero.

	ALLERRORS[cycle][0] = cycle;
	ALLERRORS[cycle][1] = n_active_cells;
	ALLERRORS[cycle][2] = n_dofs;
	//	ALLERRORS[cycle][3] = cell_diameter; // This is in class assemble_system
	// First vector solution is the numerical solution; second is the exact
	// I believe this is already all set, except for the exact solution.

	Vector<float> difference_per_cell (triangulation_new.n_active_cells());
	VectorTools::integrate_difference (dof_handler_new,
			solution_new,
			ExactSolution<dim>(),
			difference_per_cell,
			QGauss<dim>(3),
			VectorTools::L2_norm
			,&weight_values
	);
	const double L2_error = difference_per_cell.l2_norm();
	ALLERRORS[cycle][4] = L2_error;
	VectorTools::integrate_difference (dof_handler_new,
			solution_new,
			ExactSolution<dim>(),
			difference_per_cell,
			QGauss<dim>(3),
			VectorTools::H1_seminorm
			,&weight_values
	);
	const double H1_error = difference_per_cell.l2_norm();
	//    H1_norm :The square of this norm is the square of the L2_norm plus the square of the H1_seminorm.


	ALLERRORS[cycle][5] = H1_error;

	const QTrapez<1>     q_trapez;
	const QIterated<dim> q_iterated (q_trapez, 5);
	VectorTools::integrate_difference (dof_handler_new,
			solution_new,
			ExactSolution<dim>(),
			difference_per_cell,
			q_iterated,
			VectorTools::Linfty_norm
			,&weight_values
	);
	const double Linfty_error = difference_per_cell.linfty_norm();
	ALLERRORS[cycle][6] = Linfty_error;

	cout << "******** Errors: ******** \n" ;
	cout << "L2 error: " << L2_error << endl;
	cout << "H1 error: " << H1_error << endl;
	cout << "Linfty error: " << Linfty_error << endl;

	{
		// 3D Output of DIFFERENCE SOLUTION (per dof), one for each cycle of global refinement
		for (unsigned int i=0; i<solution_new.size(); ++i)
		{
			difference_solution[i] = exact_solution[i] - solution_new[i];
		}
		DataOut<dim> data_out_difference;
		data_out_difference.attach_dof_handler (dof_handler_new);
		data_out_difference.add_data_vector (difference_solution, "difference_solution");
		data_out_difference.build_patches ();

		std::string filename_exact = save_to_folder + "/difference_solution-ref-";
		filename_exact += ('0' + total_refinement);
		filename_exact += ".vtk";
		std::ofstream output_exact(filename_exact.c_str());
		data_out_difference.write_vtk(output_exact);
	}
	{
		// 3D Output of DIFFERENCE SOLUTION (cell wise), one for each cycle
		// of global refinement//
		DataOut<dim> data_out_difference;
		data_out_difference.attach_dof_handler (dof_handler_new);
		data_out_difference.add_data_vector (difference_per_cell, "difference_solution_per_cell");
		data_out_difference.build_patches ();
//		std::string filename_difference = "difference_solution_per_cell-";
		std::string filename_difference = save_to_folder + "/difference_solution_per_cell-ref-";
//		filename_difference += ('0' + cycle);
		filename_difference += ('0' + total_refinement);
		filename_difference += ".vtk";
		std::ofstream output_difference(filename_difference.c_str());
		data_out_difference.write_vtk(output_difference);
	}
}
template <int dim>
void PoissonProblem<dim>::output_results_interpolated ()
{
	//	Output CutTriangulation as VTK file.
	Write_VTK Obj_Write_VTK(CutTriangulation);
	Obj_Write_VTK.OrganizeVertices();

	{
		std::string filename_new = save_to_folder + "/interpolated_solution-ref-";
		filename_new += ('0' + total_refinement);
		filename_new += ".vtk";
		std::ofstream output(filename_new.c_str());
		Obj_Write_VTK.InterpolateSolution(dof_handler_new,solution_new);
		std::string solution_name =  "interpolated_ubulk";
		Obj_Write_VTK.WriteFile(output,solution_name);
	}
	{
		std::string filename_new = save_to_folder + "/exact_interpolated_solution-ref-";
		filename_new += ('0' + total_refinement);
		filename_new += ".vtk";
		std::ofstream output(filename_new.c_str());
		Obj_Write_VTK.InterpolateSolution(dof_handler_new,exact_solution);
		std::string solution_name =  "exact_interpolated_ubulk";
		Obj_Write_VTK.WriteFile(output,solution_name);
	}
}

template <int dim>
void PoissonProblem<dim>::run ()
{
	const clock_t begin_time = clock();
	cycle = 0;
	n_cycles = 2;
	save_to_folder = "for_publishing"; // Multiplied constants by cell size.
	ALLERRORS.reinit(n_cycles,8);

//	while (cycle < n_cycles) // In case one wants to run several cycles of refinement
	{
		make_grid (); // This will make the first grid (cycle = 0) and then just refine.
		initialize_levelset();
		output_results_levelset();
		get_new_triangulation ();

		dof_handler_new.initialize (triangulation_new,fe);

		setup_system();
		initialize_levelset_new();
		output_results_levelset_new();
		create_new_cut_cell_mesh();
		assemble_system();
		solve ();
		output_results ();
		process_solution();
		output_results_interpolated();
//		cycle++;
	}

	// Write all error info. results in a .txt file (after all cycles are completed)
	//	Order of output:
	//	Cycle|n_active_cells|n_dofs|cell_diameter|L2 Error|H1 Error|Linf Error|L1 Cond.Number

	std::ofstream all_solutions;

//	all_solutions.open("ErrorEvaluation.txt");
	for(int unsigned i = 0; i < n_cycles; ++i) {
		for(int unsigned j = 0; j < ALLERRORS.size()[1]; ++j) {
			all_solutions << ALLERRORS[i][j] << ',';
		}
		all_solutions << std::endl;
	}
	all_solutions.close();

//	Output only (3) log (cell_diameter) | (4) log(L2 Error)|(5)log(H1 Error) |(6) log(Linf Error)
//	for easy plotting in Gnuplot.
//	std::ofstream all_solutions_new;
//	all_solutions_new.open("ErrorEvaluationGnuplot.txt");
//	for(int unsigned i = 0; i < n_cycles; ++i) {
//			for(int unsigned j = 3; j < 7; ++j) {
//				all_solutions_new << log(ALLERRORS[i][j]) << ' ';
//			}
//			all_solutions_new << std::endl;
//		}
//	all_solutions_new.close();

	std::cout << "END of RUN: Successful \n";
	 cout << "Time elapsed: \n";
	 std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;
	 cout << endl;
}
} // End namespace

int main ()
{
	using namespace dealii;
	using namespace cut_cell_method;
	const unsigned int dim = 3;

	PoissonProblem<dim> Obj_PoissonProblem;
	Obj_PoissonProblem.run ();



	return 0;
}
