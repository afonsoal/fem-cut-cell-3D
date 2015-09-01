/*
 * CutCell_Integration_3D.h
 *
 *  Created on: Aug 12, 2015
 *      Author: afonsoal
 */

#ifndef CUTCELL_INTEGRATION_3D_H_
#define CUTCELL_INTEGRATION_3D_H_
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/NewCell_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/NewMesh_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/NewFace_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/Write_VTK.h"

#include <deal.II/dofs/dof_renumbering.h>
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/Testing_3D/Polynomials3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/Testing_3D/polymul.h"

using namespace dealii;
class CutCell_Integration_3D {
public:
	CutCell_Integration_3D(std::vector<NewCell_3D> const & CutTriangulation, const int dofs_per_cell);
	virtual ~CutCell_Integration_3D();
	void initialize_cell_variables(const int cell_it_, const double jacobian_determinant_, FullMatrix<double> const & jacobian_inverse_);
	void initialize_face_variables (const int face_it);
	void initialize_line_variables (const int line_it_);
	void compute_polynomials();
//	void compute_bilinear_polynomial();
	FullMatrix<double> compute_bilinear_polynomial();
	FullMatrix<double> compute_terms_c_d(const double gamma_D, const double cell_diameter);
	Vector<double> compute_rhs(const double f_B);

	double compute_bilinear_polynomial(const int dof_i, const int dof_j);
	double compute_terms_c_d(const int dof_i, const int dof_j, const double gamma_D, const double cell_diameter);
	double compute_rhs(const int dof_i, const double f_B);
	double compute_beltrami_bilinear(const int dof_i, const int dof_j);

	void compute_beltrami_polynomials();
	void compute_beltrami_polynomials_2();
	void compute_beltrami_polynomials_3();

	double compute_coupling_term(const int dof_i, const int dof_j,
			const int multiplier_us_i,
			const int multiplier_us_j,
			const int multiplier_ub_i,
			const int multiplier_ub_j, const double b_B, const double b_S);

	double compute_surface_rhs_x2(const int dof_i, const double b_S);
	double compute_surface_rhs_12xyz(const int dof_i, const double b_S);
	double compute_surface_rhs(const int dof_i, const double b_S);
	double compute_surface_rhs_2(const int dof_i, const double b_S);
	double compute_surface_rhs_3(const int dof_i, const double b_S);
	double compute_surface_rhs_4(const int dof_i, const double b_S);
	double constraint_vector(const int dof_i);

	void initialize_cell_variables_exc_ji(const int cell_it_, const double jacobian_determinant_);
	void initialize_ji(FullMatrix<double> const & jacobian_inverse_);
private:
	enum { X = 0, Y = 1, Z = 2};
	static const int dim = 3, degree = 3;

	void initialize_quadrature();

	std::vector<NewCell_3D> CutTriangulation;
	int dofs_per_cell;
	std::vector<std::vector <double> > Matrix_Coefficients;
	std::vector<NPPolynomials::Polynomials3D<dim,dim> >vector_poly3D;

	std::vector<double> line_quadrature_weights;
	std::vector<Point<dim-2> > line_quadrature_points;
	unsigned int   line_n_q_points;

//	typedef std::vector< polynomial< double,3,degree*2 > >  vector_poly;
	//	std::vector< vector_poly > bilinear_poly;
			std::vector<std::vector< polynomial< double,3,degree*2 > > > bilinear_poly;
			std::vector<std::vector< polynomial< double,3,degree*2 > > > beltrami_bilinear_poly;
//																														old: p1_F_vector
	std::vector<std::vector<std::vector<polynomial<double,3,degree*2+1>>>> bilinear_F_vector;
	std::vector<std::vector<std::vector<polynomial<double,3,degree*2+1>>>> beltrami_bilinear_F_vector;


	std::vector<std::vector<polynomial<double,3,degree> >> /*const & */grad_phi_i;
//	std::vector<std::vector<polynomial<double,3,degree> >>/*const & */grad_phi_j;
	std::vector< std::vector<polynomial<double,3,degree> > > J_grad_phi_i;
//	std::vector< std::vector<polynomial<double,3,degree> > > J_grad_phi_j;
	std::vector<std::vector<polynomial<double,3,degree+1>>> phi_i_F_vector;


	FullMatrix<double> jacobian_inverse;
	double w, jacobian_determinant, line_length_projection;
	int gamma_, alfa, beta;
	int cell_it, face_it, line_it;
	Point <dim> face_normal;

	Point <dim-1> line_norm_projection;

};

#endif /* CUTCELL_INTEGRATION_3D_H_ */


