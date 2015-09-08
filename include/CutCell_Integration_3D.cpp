/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * CutCell_Integration_3D.cpp
 *
 * This file is part of the fem-cut-cell-3D software, built using the deal.ii
 * library. You are free to use it under the MIT License as described below.
 *
 *
 * The MIT License (MIT)
 * Copyright (c) <2015> <Afonso Alborghetti Londero>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.*/
/* ---------------------------------------------------------------------
 * deal.ii license information
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 */

#include "CutCell_Integration_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/NewCell_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/NewMesh_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/NewFace_3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/Write_VTK.h"

#include <deal.II/dofs/dof_renumbering.h>
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/Testing_3D/Polynomials3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/Testing_3D/polymul.h"
#include <deal.II/base/quadrature_lib.h>

using namespace dealii;
using std::cout;
using std::endl;
CutCell_Integration_3D::~CutCell_Integration_3D() {
	// TODO Auto-generated destructor stub
}

CutCell_Integration_3D::CutCell_Integration_3D (std::vector<NewCell_3D> const & CutTriangulation, const int dofs_per_cell)
	: CutTriangulation(CutTriangulation), dofs_per_cell(dofs_per_cell)
{
	initialize_quadrature();
//						    c0       c1    c2   c3              c4    c5            c6             								       c7
		/*std::vector<std::vector <double> > */
	Matrix_Coefficients =
				{ /*0*/ {1,   -1.0,   -1,   -1,   0.0 ,    1,    1,     0.0,    1 ,      0.0,0.0,0.0,0.0,0.0,		-1,  0.0,0.0,0.0,0.0,0.0},
				/*1*/   {0,    1.0,    0,    0,   0.0 ,   -1.0, -1.0, 0.0,    0 ,     0.0,0.0,0.0,0.0,0.0,		1,   0.0,0.0,0.0,0.0,0.0},
				/*2 */  {0,    0.0,    1,    0,   0.0 ,   -1,    0,     0.0,     -1 ,     0.0,0.0,0.0,0.0,0.0,		1,    0.0,0.0,0.0,0.0,0.0},
				/*3*/   {0,    0.0,    0,    0,   0.0 ,    1,    0,     0.0,     0 ,     0.0,0.0,0.0,0.0,0.0,		-1,  0.0,0.0,0.0,0.0,0.0},
				/*4*/   {0,    0.0,    0,    1,   0.0 ,    0,   -1,     0.0,    -1 ,     0.0,0.0,0.0,0.0,0.0,		1,   0.0,0.0,0.0,0.0,0.0},
				/*5*/   {0,    0.0,    0,    0,   0.0 ,    0,    1,     0.0,      0,      0.0,0.0,0.0,0.0,0.0,		-1,   0.0,0.0,0.0,0.0,0.0},
				/*6*/   {0,    0.0,    0,    0,   0.0 ,    0,    0,     0.0,      1,      0.0,0.0,0.0,0.0,0.0,		-1,   0.0,0.0,0.0,0.0,0.0},
				/*7*/   {0,    0.0,    0,    0,   0.0 ,    0,    0,     0.0,      0,      0.0,0.0,0.0,0.0,0.0,		1,   0.0,0.0,0.0,0.0,0.0}};
//		std::vector<NPPolynomials::Polynomials3D<dim,dim> >vector_poly3D;
		vector_poly3D.reserve(Matrix_Coefficients.size());
		for (unsigned int i = 0; i< Matrix_Coefficients.size(); ++i)
		{
			NPPolynomials::Polynomials3D<dim,dim> poly3D(Matrix_Coefficients[i]);
			vector_poly3D.push_back(poly3D);
		}
		// needs jacobian_inverse!
//		Therefore, needs to be computed for EVERY cell! Meaning that what this implementation will not save any time, except if
//		one considers constant jacobian_inverse (as it is in uniform meshes).
//		compute_polynomials();
}
void CutCell_Integration_3D::initialize_quadrature()
{
	QGauss<dim-2> line_quadrature_formula(2);
	line_n_q_points = line_quadrature_formula.size();
	line_quadrature_weights = line_quadrature_formula.get_weights();
	 line_quadrature_points = line_quadrature_formula.get_points();
}

//void CutCell_Integration_3D::initialize_jacobian_inverse(FullMatrix<double> const & jacobian_inverse_)
//{
//	jacobian_inverse = jacobian_inverse_;
//}

void CutCell_Integration_3D::initialize_cell_variables(const int cell_it_, const double jacobian_determinant_,
		FullMatrix<double> const & jacobian_inverse_)
{
	cell_it = cell_it_;
	jacobian_determinant = jacobian_determinant_;
//	jacobian_determinant = CutTriangulation[cell_it].GetCellVolume();
//	jacobian_determinant = jacobian_determinant_*CutTriangulation[cell_it].GetCellVolume();
	jacobian_inverse.reinit(3,3);
	jacobian_inverse = jacobian_inverse_;
	compute_polynomials();
}
// Doesn't really make a difference, because the limiting part is compute_bilinear_polynomial (0.06 s for each line (and face) ~ 12 per cell,
// whilst compute_poly only takes 0.008 per cell))
void CutCell_Integration_3D::initialize_cell_variables_exc_ji(const int cell_it_, const double jacobian_determinant_)
{
	cell_it = cell_it_;
	jacobian_determinant = jacobian_determinant_;
}
void CutCell_Integration_3D::initialize_ji(FullMatrix<double> const & jacobian_inverse_)
{
	jacobian_inverse.reinit(3,3);
	jacobian_inverse = jacobian_inverse_;
	compute_polynomials();
}

void CutCell_Integration_3D::initialize_face_variables (const int face_it_)
{
	face_it = face_it_;
	face_normal = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].face_normal_vector;
	w = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].w;
	gamma_ = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].gamma_;
	alfa = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].alfa;
	beta = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].beta;

	assert ((fabs(face_normal[gamma_]) != 0) && "face_normal[gamma_] = 0 ");
}

void CutCell_Integration_3D::initialize_line_variables (const int line_it_)
{
	line_it = line_it_;
	line_norm_projection =
			CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].normal_projection;
	line_length_projection =
			CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].length_projection;
}

// Integrate term (grad phi_i, grad phi_j)
void CutCell_Integration_3D::compute_polynomials()
{
//	std::vector<std::vector< polynomial< double,3,degree*2 > > > bilinear_poly;
//		bilinear_poly.reserve(8);
//	p1_F_vector.reserve(8);
		bilinear_poly.resize(8);
		bilinear_F_vector.resize(8);
			for (unsigned int i = 0; i < 8; ++i)
			{
				bilinear_poly[i].reserve(8);
				bilinear_F_vector[i].resize(8);
			}

			grad_phi_i.reserve(8);
			J_grad_phi_i.reserve(8);

			for (unsigned int i = 0; i < 8; ++i)
			{
				grad_phi_i[i].reserve(3);
				J_grad_phi_i[i].reserve(3);
			}

			for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
			{
				polynomial<double,3,degree> const & phi_i = vector_poly3D[dof_i].GetPolymul();
				std::vector<polynomial<double,3,degree> > const & temp_grad_phi_i = NPPolynomials::Grad(phi_i);
				grad_phi_i.push_back(temp_grad_phi_i);
				std::vector<polynomial<double,3,degree> > const & temp_J_grad_phi_i
				= NPPolynomials::MatrixDotVectorPoly(jacobian_inverse,grad_phi_i[dof_i]);
				J_grad_phi_i.push_back(temp_J_grad_phi_i);

				// used on the computation of RHS
				std::vector<polynomial<double,3,degree+1>> p1;
				NPPolynomials::CompF_vector(phi_i, p1);
				phi_i_F_vector.push_back(p1);
			}

			for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
			{
				for (unsigned int dof_j = 0; dof_j < 8; ++dof_j)
				{
					polynomial<double,3,degree*2> p1; p1 = 0.0;
					NPPolynomials::DotMult(p1,J_grad_phi_i[dof_i],J_grad_phi_i[dof_j]);
					bilinear_poly[dof_i].push_back(p1);
					NPPolynomials::CompF_vector(p1, bilinear_F_vector[dof_i][dof_j]);
				}
			}
}
FullMatrix<double> CutCell_Integration_3D::compute_bilinear_polynomial()
{
	FullMatrix<double> cell_matrix; cell_matrix.reinit(8,8);
		for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
		{
			for (unsigned int dof_j = 0; dof_j < 8; ++dof_j)
			{
				polynomial<double,3,degree> p1; p1 = 0.0;
				p1.reinit(bilinear_poly[dof_i][dof_j]);
				p1*=jacobian_determinant;

				polynomial<double,3,degree*2+1> p1_g; p1_g = 0.0;

				for (unsigned int coord = 0; coord<dim; ++coord)
					p1_g += bilinear_F_vector[dof_i][dof_j][coord]*face_normal[coord];

				p1_g*=jacobian_determinant;

				polynomial<double,3-1,(degree*2+1)*2> p1_g_projected =
						NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

				std::vector<polynomial<double,3-1,(degree*2+1)*2+1>> p1_H_vector;
				NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);
				//////////////////////////////////////////////////////////////////////
				for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
				{
					double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
							(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
									CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
									line_quadrature_points[q_point](0);
					double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
							(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
									CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
									line_quadrature_points[q_point](0);
					Point<2> point (alfa_t,beta_t);
					Point<2> H_t;
					H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
					H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
					H_t[X] *= line_norm_projection[X]*line_length_projection; // * n_x*1/3
					H_t[Y] *= line_norm_projection[Y]*line_length_projection; // * n_x*1/3
					H_t *= (1.0/fabs(face_normal[gamma_]));
					double H_t_multiplied = H_t[X] + H_t[Y];
					H_t_multiplied *= line_quadrature_weights[q_point];
					cell_matrix(dof_i,dof_j) +=H_t_multiplied;
				} // end q points
			} // end dof i
		} // end dof j

	return cell_matrix;
}
double CutCell_Integration_3D::compute_bilinear_polynomial(const int dof_i, const int dof_j)
{
//	FullMatrix<double> cell_matrix; cell_matrix.reinit(8,8);

//	polynomial<double,3,degree> p1; p1 = 0.0;
//					p1.reinit(bilinear_poly[dof_i][dof_j]);
//					p1*=jacobian_determinant;

					polynomial<double,3,degree*2+1> p1_g; p1_g = 0.0;

					for (unsigned int coord = 0; coord<dim; ++coord)
						p1_g += bilinear_F_vector[dof_i][dof_j][coord]*face_normal[coord];

					p1_g*=jacobian_determinant;

					polynomial<double,3-1,(degree*2+1)*2> p1_g_projected =
							NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

					std::vector<polynomial<double,3-1,(degree*2+1)*2+1>> p1_H_vector;
					NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);
					//////////////////////////////////////////////////////////////////////
					double H_t_multiplied_sum = 0.0;
					for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
					{
						double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
								(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
										CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
										line_quadrature_points[q_point](0);
						double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
								(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
										CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
										line_quadrature_points[q_point](0);
						Point<2> point (alfa_t,beta_t);
						Point<2> H_t;
						H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
						H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
						H_t[X] *= line_norm_projection[X]*line_length_projection; // * n_x*1/3
						H_t[Y] *= line_norm_projection[Y]*line_length_projection; // * n_x*1/3
						H_t *= (1.0/fabs(face_normal[gamma_]));
						double H_t_multiplied = H_t[X] + H_t[Y];
						H_t_multiplied *= line_quadrature_weights[q_point];
						H_t_multiplied_sum +=H_t_multiplied;
					} // end q points
				return H_t_multiplied_sum;
//				cell_matrix(dof_i,dof_j) +=H_t_multiplied;
}

FullMatrix<double> CutCell_Integration_3D::compute_terms_c_d(const double gamma_D, const double cell_diameter)
{
	FullMatrix<double> cell_matrix; cell_matrix.reinit(8,8);
		for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
		{
			for (unsigned int dof_j = 0; dof_j < 8; ++dof_j)
			{
				polynomial<double,3,degree> const & phi_i = vector_poly3D[dof_i].GetPolymul();
				polynomial<double,3,degree> const & phi_j = vector_poly3D[dof_j].GetPolymul();

				//				std::vector<polynomial<double,3,degree> > const & grad_phi_i = NPPolynomials::Grad(phi_i);
				//				std::vector<polynomial<double,3,degree> > const & J_grad_phi_i
				//				= NPPolynomials::MatrixDotVectorPoly(jacobian_inverse,grad_phi_i);
				//
				//				std::vector<polynomial<double,3,degree> > const & grad_phi_j = NPPolynomials::Grad(phi_j);
				//				std::vector<polynomial<double,3,degree> > const & J_grad_phi_j
				//				= NPPolynomials::MatrixDotVectorPoly(jacobian_inverse,grad_phi_j);

				//								p1 = n_gamma*grad_phi_i
				polynomial<double,3,degree> p1; p1 = 0.0;
				for (unsigned int coord = 0; coord<dim; ++coord)
					p1 += J_grad_phi_i[dof_i][coord]*face_normal[coord];

				//								p1_c = (n_gamma.*grad_phi_i)*phi_j
				polynomial<double,3,degree*2> p1_c; p1_c = 0.0;
				polymul(p1_c,p1,phi_j);

				//	p2 = (n_gamma.*grad_phi_j)
				polynomial<double,3,degree> p2; p2 = 0.0;
				for (unsigned int coord = 0; coord<dim; ++coord)
					//							p2 += J_grad_phi_j[coord]*face_normal[coord];
					p2 += J_grad_phi_i[dof_j][coord]*face_normal[coord];

				//	 p2_c = (n_gamma.*grad_phi_j)*phi_i
				polynomial<double,3,degree*2> p2_c; p2_c = 0.0;
				polymul(p2_c,p2,phi_i);

				// p1_c = (n_gamma.*grad_phi_i)*phi_j + (n_gamma.*grad_phi_j)*phi_i
				p1_c += p2_c;
				polynomial<double,3,degree*2> p1_g; p1_g = 0.0;

				//	p1_g = - (n_gamma.*grad_phi_i)*phi_j - (n_gamma.*grad_phi_j)*phi_i
				//											NPPolynomials::ConstMult(p1_g,p1_c,-1.0);
				p1_g.reinit(p1_c);
				p1_g*=-1.0;
				// End term C

				polynomial<double,3,degree*2> p1_d; p1_d = 0.0;
				polymul(p1_d,phi_i,phi_j);

				p1_d*=gamma_D*1.0/cell_diameter;

				p1_g+=p1_d;

				polynomial<double,3-1,(degree*2)*2> p1_g_projected =
						NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

				std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
				NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

				for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
				{
					double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
							(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
									CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
									line_quadrature_points[q_point](0);
					double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
							(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
									CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
									line_quadrature_points[q_point](0);

					Point<2> point (alfa_t,beta_t);
					Point<2> H_t;
					H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
					H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
					H_t[X] *= line_norm_projection[X]*line_length_projection;
					H_t[Y] *= line_norm_projection[Y]*line_length_projection;
					H_t *= (1.0/fabs(face_normal[gamma_]));
					double H_t_multiplied = H_t[X] + H_t[Y];
					H_t_multiplied *= line_quadrature_weights[q_point];
					H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
					cell_matrix(dof_i,dof_j) +=H_t_multiplied;
				} //end q_points

			} // end dof_j
		} // end dof_i
	return cell_matrix;
}
double CutCell_Integration_3D::compute_terms_c_d(const int dof_i, const int dof_j,const double gamma_D, const double cell_diameter)
{
	assert(dof_i<8 && dof_j<8);
				polynomial<double,3,degree> const & phi_i = vector_poly3D[dof_i].GetPolymul();
				polynomial<double,3,degree> const & phi_j = vector_poly3D[dof_j].GetPolymul();

				//								p1 = n_gamma*grad_phi_i
				polynomial<double,3,degree> p1; p1 = 0.0;
				for (unsigned int coord = 0; coord<dim; ++coord)
					p1 += J_grad_phi_i[dof_i][coord]*face_normal[coord];

				//								p1_c = (n_gamma.*grad_phi_i)*phi_j
				polynomial<double,3,degree*2> p1_c; p1_c = 0.0;
				polymul(p1_c,p1,phi_j);

				//	p2 = (n_gamma.*grad_phi_j)
				polynomial<double,3,degree> p2; p2 = 0.0;
				for (unsigned int coord = 0; coord<dim; ++coord)
					//							p2 += J_grad_phi_j[coord]*face_normal[coord];
					p2 += J_grad_phi_i[dof_j][coord]*face_normal[coord];

				//	 p2_c = (n_gamma.*grad_phi_j)*phi_i
				polynomial<double,3,degree*2> p2_c; p2_c = 0.0;
				polymul(p2_c,p2,phi_i);

				// p1_c = (n_gamma.*grad_phi_i)*phi_j + (n_gamma.*grad_phi_j)*phi_i
				p1_c += p2_c;
				polynomial<double,3,degree*2> p1_g; p1_g = 0.0;

				//	p1_g = - (n_gamma.*grad_phi_i)*phi_j - (n_gamma.*grad_phi_j)*phi_i
				//											NPPolynomials::ConstMult(p1_g,p1_c,-1.0);
				p1_g.reinit(p1_c);
				p1_g*=-1.0;
				// End term C

				polynomial<double,3,degree*2> p1_d; p1_d = 0.0;
				polymul(p1_d,phi_i,phi_j);

				p1_d*=gamma_D*1.0/cell_diameter;

				p1_g+=p1_d;

				polynomial<double,3-1,(degree*2)*2> p1_g_projected =
						NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

				std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
				NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

				double H_t_multiplied_sum = 0.0;
				for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
				{
					double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
							(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
									CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
									line_quadrature_points[q_point](0);
					double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
							(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
									CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
									line_quadrature_points[q_point](0);

					Point<2> point (alfa_t,beta_t);
					Point<2> H_t;
					H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
					H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
					H_t[X] *= line_norm_projection[X]*line_length_projection;
					H_t[Y] *= line_norm_projection[Y]*line_length_projection;
					H_t *= (1.0/fabs(face_normal[gamma_]));
					double H_t_multiplied = H_t[X] + H_t[Y];
					H_t_multiplied *= line_quadrature_weights[q_point];
					H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
					H_t_multiplied_sum +=H_t_multiplied;
				} //end q_points

	return H_t_multiplied_sum;
}

Vector<double> CutCell_Integration_3D::compute_rhs(const double f_B)
{

//	Vector<double> cell_rhs; cell_rhs.reinit(8);
	Vector<double> cell_rhs(dofs_per_cell); cell_rhs = 0.0;
		for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
		{
//			polynomial<double,3,degree> const & phi_i = vector_poly3D[dof_i].GetPolymul();
//			std::vector<polynomial<double,3,degree+1>> phi_i_F_vector;
//			NPPolynomials::CompF_vector(phi_i, phi_i_F_vector);

			polynomial<double,3,degree+1> p1_g;
			p1_g = 0.0;

			for (unsigned int coord = 0; coord<dim; ++coord)
				p1_g += phi_i_F_vector[dof_i][coord]*face_normal[coord];

			p1_g*=jacobian_determinant;

			polynomial<double,3-1,(degree+1)*2> p1_g_projected =
					NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

			std::vector<polynomial<double,3-1,(degree+1)*2+1>> p1_H_vector;
			NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

			for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
			{
				double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
						(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
								CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
								line_quadrature_points[q_point](0);
				double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
						(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
								CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
								line_quadrature_points[q_point](0);

			Point<2> point (alfa_t,beta_t);
			Point<2> H_t;

			H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
			H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
			H_t[X] *= line_norm_projection[X]*line_length_projection; // * n_x*1/3
			H_t[Y] *= line_norm_projection[Y]*line_length_projection; // * n_x*1/3
			H_t 	*= (1.0/fabs(face_normal[gamma_]));
			double H_t_multiplied = H_t[X] + H_t[Y];
			H_t_multiplied *= line_quadrature_weights[q_point];
			cell_rhs(dof_i) += f_B*H_t_multiplied;
		} // end q_points
	} // end dof_i
	return cell_rhs;
}

double CutCell_Integration_3D::compute_rhs(const int dof_i, const double f_B)
{
		assert(dof_i<8);
			polynomial<double,3,degree+1> p1_g;
			p1_g = 0.0;

			for (unsigned int coord = 0; coord<dim; ++coord)
				p1_g += phi_i_F_vector[dof_i][coord]*face_normal[coord];

			p1_g*=jacobian_determinant;

			polynomial<double,3-1,(degree+1)*2> p1_g_projected =
					NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

			std::vector<polynomial<double,3-1,(degree+1)*2+1>> p1_H_vector;
			NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

			double H_t_multiplied_sum = 0.0;
			for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
			{
				double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
						(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
								CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
								line_quadrature_points[q_point](0);
				double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
						(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
								CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
								line_quadrature_points[q_point](0);

			Point<2> point (alfa_t,beta_t);
			Point<2> H_t;

			H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
			H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
			H_t[X] *= line_norm_projection[X]*line_length_projection; // * n_x*1/3
			H_t[Y] *= line_norm_projection[Y]*line_length_projection; // * n_x*1/3
			H_t 	*= (1.0/fabs(face_normal[gamma_]));
			double H_t_multiplied = H_t[X] + H_t[Y];
			H_t_multiplied *= line_quadrature_weights[q_point];
			H_t_multiplied_sum += f_B*H_t_multiplied;
		} // end q_points
	return H_t_multiplied_sum;
}

// Has to be called for every face because it needs face_normal.
void CutCell_Integration_3D::compute_beltrami_polynomials()
{
	beltrami_bilinear_poly.clear();
	beltrami_bilinear_poly.resize(dofs_per_cell);

	FullMatrix<double> beltrami_operator_0;
	beltrami_operator_0.reinit(dim,dim);
	beltrami_operator_0 = IdentityMatrix(dim);
	Tensor<2,dim> n_outer_n;
	outer_product(n_outer_n, face_normal,face_normal);
	FullMatrix<double> matrix_n_outer_n(dim,dim);
	matrix_n_outer_n.copy_from(n_outer_n);
	FullMatrix<double> beltrami_operator(dim,dim);
	beltrami_operator.copy_from(beltrami_operator_0);
//	beltrami_operator.add(1.0,matrix_n_outer_n);
	beltrami_operator.add(-1.0,matrix_n_outer_n);

	for (unsigned int dof_i = 0; dof_i < dofs_per_cell; ++dof_i)
		beltrami_bilinear_poly[dof_i].reserve(dofs_per_cell);

	for (unsigned int dof_i = 0; dof_i < dofs_per_cell; ++dof_i)
	{
		std::vector<polynomial<double,dim,degree>>  const & temp0 =
					NPPolynomials::MatrixDotVectorPoly(beltrami_operator,J_grad_phi_i[dof_i]);
		for (unsigned int dof_j = 0; dof_j < dofs_per_cell; ++dof_j)
		{
			 std::vector<polynomial<double,dim,degree>>  const & temp1 =
			 			NPPolynomials::MatrixDotVectorPoly(beltrami_operator,J_grad_phi_i[dof_j]);

			 polynomial<double,3,degree*2> p1; p1 = 0.0;
			NPPolynomials::DotMult(p1,temp0,temp1);
			beltrami_bilinear_poly[dof_i].push_back(p1);
//			NPPolynomials::CompF_vector(p1, beltrami_bilinear_F_vector[dof_i][dof_j]);
		}
		assert(beltrami_bilinear_poly[dof_i].size() == dofs_per_cell);
	}

}
// Multiply J_inverse at the end. Yields same result.
void CutCell_Integration_3D::compute_beltrami_polynomials_2()
{
	beltrami_bilinear_poly.clear();
	beltrami_bilinear_poly.resize(dofs_per_cell);

	FullMatrix<double> beltrami_operator_0;
	beltrami_operator_0.reinit(dim,dim);
	beltrami_operator_0 = IdentityMatrix(dim);
	Tensor<2,dim> n_outer_n;
	outer_product(n_outer_n, face_normal,face_normal);
	FullMatrix<double> matrix_n_outer_n(dim,dim);
	matrix_n_outer_n.copy_from(n_outer_n);
	FullMatrix<double> beltrami_operator(dim,dim);
	beltrami_operator.copy_from(beltrami_operator_0);
//	beltrami_operator.add(1.0,matrix_n_outer_n);
	beltrami_operator.add(-1.0,matrix_n_outer_n);
//	beltrami_operator.add(0.0,matrix_n_outer_n);

	for (unsigned int dof_i = 0; dof_i < dofs_per_cell; ++dof_i)
		beltrami_bilinear_poly[dof_i].reserve(dofs_per_cell);

	for (unsigned int dof_i = 0; dof_i < dofs_per_cell; ++dof_i)
	{
		std::vector<polynomial<double,dim,degree>>  temp0 =
					NPPolynomials::MatrixDotVectorPoly(beltrami_operator,/*J_*/grad_phi_i[dof_i]);
		temp0 = NPPolynomials::MatrixDotVectorPoly(jacobian_inverse,temp0);

		for (unsigned int dof_j = 0; dof_j < dofs_per_cell; ++dof_j)
		{
			 std::vector<polynomial<double,dim,degree>>  temp1 =
			 			NPPolynomials::MatrixDotVectorPoly(beltrami_operator,/*J_*/grad_phi_i[dof_j]);
			 temp1 = NPPolynomials::MatrixDotVectorPoly(jacobian_inverse,temp1);


			 polynomial<double,3,degree*2> p1; p1 = 0.0;
			NPPolynomials::DotMult(p1,temp0,temp1);
			beltrami_bilinear_poly[dof_i].push_back(p1);
//			NPPolynomials::CompF_vector(p1, beltrami_bilinear_F_vector[dof_i][dof_j]);
		}
		assert(beltrami_bilinear_poly[dof_i].size() == dofs_per_cell);
	}

}

// Alternatively, using approach by https://www.igpm.rwth-aachen.de/Download/reports/pdf/IGPM342_k.pdf
// ... which is of course the same thing as using identity matrix, etc.
//Yields same result as compute_beltrami_polynomials.
void CutCell_Integration_3D::compute_beltrami_polynomials_3()
{
	beltrami_bilinear_poly.clear();
	beltrami_bilinear_poly.resize(8);

	for (unsigned int i = 0; i < 8; ++i)
		beltrami_bilinear_poly[i].reserve(8);

	std::vector<std::vector<polynomial<double,3,degree> >> beltrami_grad_phi_i;
	std::vector< std::vector<polynomial<double,3,degree> > > J_beltrami_grad_phi_i;

	beltrami_grad_phi_i.resize(8);
	J_beltrami_grad_phi_i.resize(8);

			for (unsigned int i = 0; i < 8; ++i)
			{
				beltrami_grad_phi_i[i].reserve(3);
				J_beltrami_grad_phi_i[i].reserve(3);
			}

			for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
			{
				polynomial<double,3,degree> p0; p0 = 0.0;
				for (unsigned int coord = 0; coord<dim; ++coord)
					p0 += J_grad_phi_i[dof_i][coord]*face_normal[coord];

				p0*=-1.;
				std::vector<polynomial<double,3,degree>> p1;
				for (unsigned int coord = 0; coord<dim; ++coord)
				{
					polynomial<double,3,degree> temp; temp = 0.;
					temp.reinit(p0); temp*=face_normal[coord];
					p1.push_back(temp);
				}

				std::vector<polynomial<double,3,degree>> p2/*(dim)*/;
	//			p2.resize(dim);
				for (unsigned int coord = 0; coord<dim; ++coord)
				{
					polynomial<double,3,degree> temp; temp = 0.;
					temp.reinit(J_grad_phi_i[dof_i][coord]);
					temp+=p1[coord];
					J_beltrami_grad_phi_i[dof_i].push_back(temp);
				}
			}

			for (unsigned int dof_i = 0; dof_i < 8; ++dof_i)
			{
				for (unsigned int dof_j = 0; dof_j < 8; ++dof_j)
				{
					polynomial<double,3,degree*2> p1; p1 = 0.0;
					NPPolynomials::DotMult(p1,J_beltrami_grad_phi_i[dof_i],J_beltrami_grad_phi_i[dof_j]);
					beltrami_bilinear_poly[dof_i].push_back(p1);
//					NPPolynomials::CompF_vector(p1, bilinear_F_vector[dof_i][dof_j]);
				}
			}
}

double CutCell_Integration_3D::compute_beltrami_bilinear(const int dof_i, const int dof_j)
{


//	for (unsigned int coord = 0; coord<dim; ++coord)
//		p1_g +=  beltrami_bilinear_F_vector[dof_i][dof_j][coord]*face_normal[coord];

	polynomial<double,3,degree*2> p1_g; p1_g = 0.0;
	p1_g.reinit(beltrami_bilinear_poly[dof_i][dof_j]);

	polynomial<double,3-1,(degree*2)*2> p1_g_projected =
			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

	std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

	double H_t_multiplied_sum = 0.0;
	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
	{
		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
						line_quadrature_points[q_point](0);
		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
						line_quadrature_points[q_point](0);

		Point<2> point (alfa_t,beta_t);
		Point<2> H_t;
		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
		H_t[X] *= line_norm_projection[X]*line_length_projection;
		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
		H_t *= (1.0/fabs(face_normal[gamma_]));
		double H_t_multiplied = H_t[X] + H_t[Y];
		H_t_multiplied *= line_quadrature_weights[q_point];
		H_t_multiplied*=pow(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea(),0.5);
		H_t_multiplied_sum +=H_t_multiplied;
	} //end q_points
	return H_t_multiplied_sum;
}

double CutCell_Integration_3D::compute_coupling_term(const int dof_i, const int dof_j,
		const int multiplier_us_i,
		const int multiplier_us_j,
		const int multiplier_ub_i,
		const int multiplier_ub_j, const double b_B, const double b_S)
{

	polynomial<double,3,degree> const & phi_i = vector_poly3D[dof_i].GetPolymul();
	polynomial<double,3,degree> const & phi_j = vector_poly3D[dof_j].GetPolymul();

	polynomial<double,dim,degree*2> phi_ij; phi_ij = 0.0;
	polymul(phi_ij,phi_i,phi_j);

	polynomial<double,dim,degree*2> p0; p0 = 0.0;
	polynomial<double,dim,degree*2> p1; p1 = 0.0;
	polynomial<double,dim,degree*2> p2; p2 = 0.0;
	polynomial<double,dim,degree*2> p3; p3 = 0.0;

	p0.reinit(phi_ij);
	p0*=b_B*b_B*multiplier_ub_i* multiplier_ub_j;

	p1.reinit(phi_ij);
	p1*=-b_B*b_S*multiplier_ub_i* multiplier_us_j;

	p2.reinit(phi_ij);
	p2*=-b_B*b_S*multiplier_us_i* multiplier_ub_j;

	p3.reinit(phi_ij);
	p3*=b_B*b_S*multiplier_us_i* multiplier_us_j;

	polynomial<double,3,degree*2> p1_g; p1_g = 0.0;

	p1_g+=p0+p1+p2+p3;

		polynomial<double,3-1,(degree*2)*2> p1_g_projected =
				NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

		std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
		NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

		double H_t_multiplied_sum = 0.0;
		for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
		{
			double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
					(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
							CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
							line_quadrature_points[q_point](0);
			double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
					(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
							CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
							line_quadrature_points[q_point](0);

			Point<2> point (alfa_t,beta_t);
			Point<2> H_t;
			H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
			H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
			H_t[X] *= line_norm_projection[X]*line_length_projection;
			H_t[Y] *= line_norm_projection[Y]*line_length_projection;
			H_t *= (1.0/fabs(face_normal[gamma_]));
			double H_t_multiplied = H_t[X] + H_t[Y];
			H_t_multiplied *= line_quadrature_weights[q_point];
			H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
			H_t_multiplied_sum +=H_t_multiplied;
		} //end q_points
		return H_t_multiplied_sum;
}

//Compute surface RHS : f_S*phi_i, where f_S=2 (...) s.t exact_solution = x^2+y^2+z^2.
double CutCell_Integration_3D::compute_surface_rhs_x2(const int dof_i, const double b_S)
{

	polynomial<double,3,degree> p1_g = vector_poly3D[dof_i].GetPolymul();

	p1_g*=2.*b_S;

	polynomial<double,3-1,(degree)*2> p1_g_projected =
			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

	std::vector<polynomial<double,3-1,(degree)*2+1>> p1_H_vector;
	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

	double H_t_multiplied_sum = 0.0;
	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
	{
		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
						line_quadrature_points[q_point](0);
		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
						line_quadrature_points[q_point](0);

		Point<2> point (alfa_t,beta_t);
		Point<2> H_t;
		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
		H_t[X] *= line_norm_projection[X]*line_length_projection;
		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
		H_t *= (1.0/fabs(face_normal[gamma_]));
		double H_t_multiplied = H_t[X] + H_t[Y];
		H_t_multiplied *= line_quadrature_weights[q_point];
		H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
		H_t_multiplied_sum +=H_t_multiplied;
	} //end q_points
	return H_t_multiplied_sum;
}
double CutCell_Integration_3D::compute_surface_rhs_12xyz(const int dof_i, const double b_S)
{


	//		Option 1. f = xyz
			polynomial<double,3,degree> f_S; f_S = 0.0;
			int expon[dim] = {1,1,1};
			f_S[f_S.term_index(expon)] = 12.;
			assert(f_S[14] == 12.0);

		polynomial<double,3,degree> phi_i = vector_poly3D[dof_i].GetPolymul();

		polynomial<double,3,degree*2> p1_g; p1_g = 0.;
		polymul(p1_g,f_S,phi_i);
		p1_g*=b_S;

	polynomial<double,3-1,(degree*2)*2> p1_g_projected =
			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

	std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

	double H_t_multiplied_sum = 0.0;
	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
	{
		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
						line_quadrature_points[q_point](0);
		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
						line_quadrature_points[q_point](0);

		Point<2> point (alfa_t,beta_t);
		Point<2> H_t;
		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
		H_t[X] *= line_norm_projection[X]*line_length_projection;
		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
		H_t *= (1.0/fabs(face_normal[gamma_]));
		double H_t_multiplied = H_t[X] + H_t[Y];
		H_t_multiplied *= line_quadrature_weights[q_point];
		H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
		H_t_multiplied_sum +=H_t_multiplied;
	} //end q_points
	return H_t_multiplied_sum;
}

//Compute surface RHS : f_S*phi_i, where f_S=12 xyz s.t exact_solution = xyz.
double CutCell_Integration_3D::compute_surface_rhs(const int dof_i, const double b_S)
{

//	Option 1. f = xyz
//	polynomial<double,3,degree> f_S; f_S = 0.0;
//	int expon[dim] = {1,1,1};
//	f_S[f_S.term_index(expon)] = 12.;
//	assert(f_S[14] == 12.0);
	//	Option 2. f = xy
	//	int expon[dim] = {1,1,0};
//	f_S[f_S.term_index(expon)] = 12.;

	//Option 3. f = 6(x+y+z), u = x^3+y^3+z^3
//	polynomial<double,3,degree> f_S; f_S = 0.0;
//		int expon[dim] = {1,0,0};
//		f_S[f_S.term_index(expon)] = -6.;
////		assert(f_S[1]==6.);
//		int expon2[dim] = {0,1,0};
//		f_S[f_S.term_index(expon2)] = -6.;
////		assert(f_S[2]==6.);
//		int expon3[dim] = {0,0,1};
//		f_S[f_S.term_index(expon3)] = -6.;
//		assert(f_S[3]==6.);

//		Option 4. f is calculated using L-B operator. u = x^3+y^3+z^3
		polynomial<double,3,degree> f_S; f_S = 0.0;
			int expon[dim] = {1,0,0};
			f_S[f_S.term_index(expon)] = -6.*pow((1-face_normal[X]*face_normal[X]),2);
			int expon2[dim] = {0,1,0};
			f_S[f_S.term_index(expon2)] = -6.*pow((1-face_normal[Y]*face_normal[Y]),2);
			int expon3[dim] = {0,0,1};
			f_S[f_S.term_index(expon3)] = -6.*pow((1-face_normal[Z]*face_normal[Z]),2);


//	assert(f_S[14] == 1.0);
	polynomial<double,3,degree*2> p1_g; p1_g = 0.0;
	polynomial<double,3,degree> phi_i = vector_poly3D[dof_i].GetPolymul();
//	polynomial<double,3,degree> p1_g; p1_g = 0.0;
	polymul(p1_g,f_S,phi_i);
//	p1_g.reinit(phi_i);

	p1_g*=b_S;

	polynomial<double,3-1,(degree*2)*2> p1_g_projected =
			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

	std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

	double H_t_multiplied_sum = 0.0;
	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
	{
		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
						line_quadrature_points[q_point](0);
		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
						line_quadrature_points[q_point](0);

		Point<2> point (alfa_t,beta_t);
		Point<2> H_t;
		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
		H_t[X] *= line_norm_projection[X]*line_length_projection;
		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
		H_t *= (1.0/fabs(face_normal[gamma_]));
		double H_t_multiplied = H_t[X] + H_t[Y];
		H_t_multiplied *= line_quadrature_weights[q_point];
		H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
		H_t_multiplied_sum +=H_t_multiplied;
	} //end q_points
	return H_t_multiplied_sum;
}

//Compute surface RHS : f_S is "constant", does not enter in the integral (doesn't work nor make sense!)
double CutCell_Integration_3D::compute_surface_rhs_4(const int dof_i, const double b_S)
{

	polynomial<double,3,degree> p1_g = vector_poly3D[dof_i].GetPolymul();

	p1_g*=b_S;

	polynomial<double,3-1,(degree)*2> p1_g_projected =
			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

	std::vector<polynomial<double,3-1,(degree)*2+1>> p1_H_vector;
	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

	double H_t_multiplied_sum = 0.0;
	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
	{
		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
						line_quadrature_points[q_point](0);
		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
						line_quadrature_points[q_point](0);

		Point<2> point (alfa_t,beta_t);
		Point<2> H_t;
		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
		H_t[X] *= line_norm_projection[X]*line_length_projection;
		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
		H_t *= (1.0/fabs(face_normal[gamma_]));
		double H_t_multiplied = H_t[X] + H_t[Y];
		H_t_multiplied *= line_quadrature_weights[q_point];
		H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
		H_t_multiplied_sum +=H_t_multiplied;
	} //end q_points
	return H_t_multiplied_sum;
}

// Use automated polynomial approach, considering J⁻¹ div  (J⁻¹ grad u) = f_S, u = x³+y³+z³
double CutCell_Integration_3D::compute_surface_rhs_2(const int dof_i, const double b_S){}
//{
//
//		polynomial<double,3,degree> u; u = 0.0;
//	//Option 4. f is calculated using L-B operator. u = x^3+y^3+z^3
//			int expon[dim] = {3,0,0};
//			u[u.term_index(expon)] = 1.;
//			int expon2[dim] = {0,3,0};
//			u[u.term_index(expon2)] = 1.;
//			int expon3[dim] = {0,0,3};
//			u[u.term_index(expon3)] = 1.;
//
//			std::vector<polynomial<double,3,degree> > const & grad_u = NPPolynomials::Grad(u);
//
//			Tensor<2,dim> n_outer_n;
//			outer_product(n_outer_n, face_normal,face_normal);
//			FullMatrix<double> matrix_n_outer_n(dim,dim);
//			matrix_n_outer_n.copy_from(n_outer_n);
//			matrix_n_outer_n*=-1.;
//
//			std::vector<polynomial<double,3,degree> > const & grad_u_nnt
//						= NPPolynomials::MatrixDotVectorPoly(matrix_n_outer_n,grad_u);
//
//			std::vector<polynomial<double,3,degree> > const grad_beltrami_u = grad_u+grad_u_nnt;
//
//			polynomial<double,3,degree> div_grad_beltrami_u;
//
//			std::vector<polynomial<double,3,degree>> temp_div_grad_u;
//			temp_div_grad_u[X] = 0.;
//			temp_div_grad_u[X] = NPPolynomials::Grad(div_grad_beltrami_u[X])[X];
//
//			temp_div_grad_u[Y] = 0.;
//			temp_div_grad_u[Y] = NPPolynomials::Grad(div_grad_beltrami_u[Y])[Y];
//
//			temp_div_grad_u[Z] = 0.;
//			temp_div_grad_u[Z] = NPPolynomials::Grad(div_grad_beltrami_u[Z])[Z];
//
//			std::vector<polynomial<double,3,degree> > div_grad_u_nnt
//						= NPPolynomials::MatrixDotVectorPoly(matrix_n_outer_n,grad_u);
//
//
//
////	assert(f_S[14] == 1.0);
//	polynomial<double,3,degree*2> p1_g; p1_g = 0.0;
//	polynomial<double,3,degree> phi_i = vector_poly3D[dof_i].GetPolymul();
//	polymul(p1_g,J_div_grad_u_nnt,phi_i);
////	p1_g.reinit(phi_i);
//
//	p1_g*=-1.*b_S;
//
//	polynomial<double,3-1,(degree*2)*2> p1_g_projected =
//			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);
//
//	std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
//	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);
//
//	double H_t_multiplied_sum = 0.0;
//	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
//	{
//		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
//				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
//						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
//						line_quadrature_points[q_point](0);
//		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
//				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
//						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
//						line_quadrature_points[q_point](0);
//
//		Point<2> point (alfa_t,beta_t);
//		Point<2> H_t;
//		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
//		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
//		H_t[X] *= line_norm_projection[X]*line_length_projection;
//		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
//		H_t *= (1.0/fabs(face_normal[gamma_]));
//		double H_t_multiplied = H_t[X] + H_t[Y];
//		H_t_multiplied *= line_quadrature_weights[q_point];
//		H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
//		H_t_multiplied_sum +=H_t_multiplied;
//	} //end q_points
//	return H_t_multiplied_sum;
//}

double CutCell_Integration_3D::compute_surface_rhs_3(const int dof_i, const double b_S)
{
	FullMatrix<double> beltrami_operator_0;
	beltrami_operator_0.reinit(dim,dim);
	beltrami_operator_0 = IdentityMatrix(dim);
	Tensor<2,dim> n_outer_n;
	outer_product(n_outer_n, face_normal,face_normal);
	FullMatrix<double> matrix_n_outer_n(dim,dim);
	matrix_n_outer_n.copy_from(n_outer_n);
	FullMatrix<double> beltrami_operator(dim,dim);
	beltrami_operator.copy_from(beltrami_operator_0);
//	beltrami_operator.add(1.0,matrix_n_outer_n);
	beltrami_operator.add(-1.0,matrix_n_outer_n);


	polynomial<double,3,degree> u; u = 0.0;
//Option 4. f is calculated using L-B operator. u = x^3+y^3+z^3
		int expon[dim] = {3,0,0};
		u[u.term_index(expon)] = -1.;
//			assert(f_S[1]==6.);
		int expon2[dim] = {0,3,0};
		u[u.term_index(expon2)] = -1.;
//			assert(f_S[2]==6.);
		int expon3[dim] = {0,0,1};
		u[u.term_index(expon3)] = -1.;
//			assert(f_S[3]==6.);

	std::vector<polynomial<double,3,degree> > const & grad_u = NPPolynomials::Grad(u);
	std::vector<polynomial<double,3,degree> > const & J_grad_u
	= NPPolynomials::MatrixDotVectorPoly(jacobian_inverse,grad_u);

	 std::vector<polynomial<double,dim,degree>>  const & J_grad_S_u_nnt =
	NPPolynomials::MatrixDotVectorPoly(beltrami_operator,J_grad_u);
		polynomial<double,3,degree> J_div_grad_beltrami_u;

		std::vector<polynomial<double,3,degree>> temp_div_J_grad_S_u(dim);
		temp_div_J_grad_S_u[X] = 0.;
		temp_div_J_grad_S_u[X] = NPPolynomials::Grad(J_grad_S_u_nnt[X])[X];

		temp_div_J_grad_S_u[Y] = 0.;
		temp_div_J_grad_S_u[Y] = NPPolynomials::Grad(J_grad_S_u_nnt[Y])[Y];

		temp_div_J_grad_S_u[Z] = 0.;
		temp_div_J_grad_S_u[Z] = NPPolynomials::Grad(J_grad_S_u_nnt[Z])[Z];



		polynomial<double,3,degree> div_J_grad_S_u = 0.;
		for (unsigned int coord = 0; coord<dim; ++coord)
			div_J_grad_S_u += temp_div_J_grad_S_u[coord];

		for (unsigned int coord = 0; coord<dim; ++coord)
			temp_div_J_grad_S_u[coord]*=-face_normal[coord]*face_normal[coord];

		for (unsigned int coord = 0; coord<dim; ++coord)
			div_J_grad_S_u += temp_div_J_grad_S_u[coord];

		if (cell_it == 0 && dof_i == 0)
        for (int i=0;i<div_J_grad_S_u.size; ++i)
        {
        	div_J_grad_S_u.exponents(i,expon);
       	 std::cout << i << " div_J_grad_S_u[i]: " << div_J_grad_S_u[i] << " expon: " << expon[0] << " " << expon[1] << " " << expon[2] << std::endl;
        }


			polynomial<double,3,degree*2> p1_g; p1_g = 0.0;
			polynomial<double,3,degree> phi_i = vector_poly3D[dof_i].GetPolymul();
			polymul(p1_g,div_J_grad_S_u,phi_i);
		//	p1_g.reinit(phi_i);

			p1_g*=b_S;
//			cout << "Test 3\n";
			polynomial<double,3-1,(degree*2)*2> p1_g_projected =
					NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

			std::vector<polynomial<double,3-1,(degree*2)*2+1>> p1_H_vector;
			NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

			double H_t_multiplied_sum = 0.0;
			for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
			{
				double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
						(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
								CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
								line_quadrature_points[q_point](0);
				double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
						(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
								CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
								line_quadrature_points[q_point](0);

				Point<2> point (alfa_t,beta_t);
				Point<2> H_t;
				H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
				H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
				H_t[X] *= line_norm_projection[X]*line_length_projection;
				H_t[Y] *= line_norm_projection[Y]*line_length_projection;
				H_t *= (1.0/fabs(face_normal[gamma_]));
				double H_t_multiplied = H_t[X] + H_t[Y];
				H_t_multiplied *= line_quadrature_weights[q_point];
				H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
				H_t_multiplied_sum +=H_t_multiplied;
			} //end q_points
			return H_t_multiplied_sum;
}

double CutCell_Integration_3D::constraint_vector(const int dof_i)
{
	polynomial<double,3,degree> p1_g; p1_g = 0.0;
	polynomial<double,3,degree> phi_i = vector_poly3D[dof_i].GetPolymul();
	p1_g.reinit(phi_i);

	polynomial<double,3-1,(degree)*2> p1_g_projected =
			NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

	std::vector<polynomial<double,3-1,(degree)*2+1>> p1_H_vector;
	NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);

	double H_t_multiplied_sum = 0.0;
	for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
	{
		double alfa_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
						line_quadrature_points[q_point](0);
		double beta_t = CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
				(CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
						CutTriangulation[cell_it].Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
						line_quadrature_points[q_point](0);

		Point<2> point (alfa_t,beta_t);
		Point<2> H_t;
		H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
		H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
		H_t[X] *= line_norm_projection[X]*line_length_projection;
		H_t[Y] *= line_norm_projection[Y]*line_length_projection;
		H_t *= (1.0/fabs(face_normal[gamma_]));
		double H_t_multiplied = H_t[X] + H_t[Y];
		H_t_multiplied *= line_quadrature_weights[q_point];
		H_t_multiplied*=CutTriangulation[cell_it].Obj_VectorNewFace[face_it].GetFaceArea();
		H_t_multiplied_sum +=H_t_multiplied;
	} //end q_points
	return H_t_multiplied_sum;
}
