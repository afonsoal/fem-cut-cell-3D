/* The fem-cut-cell-3D software
 * Author	  : Afonso Alborghetti Londero (afonsoal@gmail.com)
 * Last update:	08/Sep/2015
 *
 * NewCell_3D.cpp
 *
 * This file is part of the fem-cut-cell-3D software, built using the deal.ii
 * library. You are free to use it under the GNU Lesser General Public License
 * as described in the LICENSE File.
 * Copyright (c) <2015> <Afonso Alborghetti Londero>
 *
 * ---------------------------------------------------------------------
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

#include "NewCell_3D.h"
#include <fstream>
#include <iostream>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include "NewFace_3D.h"
#include <math.h> // fabs
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/Testing_3D/Polynomials3D.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/Testing_3D/polymul.h"

using namespace dealii;

// Constructor.
NewCell_3D::NewCell_3D(	bool _is_surface_cell)
{
	cell_volume = 0.;
	is_surface_cell = _is_surface_cell;
	if (is_surface_cell == true)
		number_of_faces = 0;

	real_faces_vector.reserve(GeometryInfo<3>::faces_per_cell);
}
void NewCell_3D::InputNewRealFaceInfo(int _face_it)
{
	real_faces_info obj_real_face;
	obj_real_face.face_it = _face_it;
	obj_real_face.is_stabilization_face = false /*_is_stabilization_face*/;
	real_faces_vector.push_back(obj_real_face);
}
// Old SetNewRealFaceStabilization
// Define the real face identified by the iterator as a face that is part of the integration of the stabilization term.
void NewCell_3D::DefStabilizationFace ( int real_face_it, bool _is_stabilization_face)
{
	real_faces_vector[real_face_it].is_stabilization_face = _is_stabilization_face;
}
// Input a new CutFace.
void NewCell_3D::InputNewFace (NewFace_3D &CutFace)
{

	CutFace.SetFaceIndex(number_of_faces);
	++number_of_faces;
	Obj_VectorNewFace.push_back(CutFace);
}

// Set the index of the cell to which this face belongs to.
void NewCell_3D::SetIndex(/*hp::*/DoFHandler<3>::active_cell_iterator _index)
{
	cell_index = _index;
}

//Compute new boundary face,
void NewCell_3D::CompBoundaryFace(NewFace_3D &CutFaceBoundary)
{
	int line_no = 0;
	bool add_1_more_line = false;
	for (int face_it = 0; face_it < number_of_faces; ++face_it) // Here the CutBoundaryFace is not counted yet.
		for (int line_it = 0; line_it <  Obj_VectorNewFace[face_it].number_of_lines/*CutFaceBoundary.number_of_lines*/; ++line_it)
	{
//			std::cout << "face_it: " << face_it << " line_it: " << line_it << " is_boundary_face: " << Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].is_boundary_face << " end" << std::endl;
			if (Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].is_boundary_line)
			{
//				Use real X0, X1.
				Point<3> X0 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0;
				Point<3> X1 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1;
				// NEW: USE mapped unit X0, X1
//				Point<3> X0 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0_unit;
//				Point<3> X1 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1_unit;
				for (unsigned int coord = 0; coord<3; ++coord)
				{
					if (fabs(X0[coord]) < pow(10,-10))
							X0[coord] = 0.0;

					if (fabs(X1[coord]) < pow(10,-10))
							X1[coord] = 0.0;
				}
				// SetCoordinates inputs newLine
				// Up to this point, I already know how many lines CutFaceBoundary has, because I sum 1 line every time I find a cut line.
				// (in the main.cpp) (why??). That's why here I do not need to add_1_more_line.
				 CutFaceBoundary.SetCoordinates(line_no,X0,X1,true /*is boundary (line)*/,add_1_more_line);
				line_no++;
			}
	}
	CutFaceBoundary.is_boundary_face = true;
//	CutFaceBoundary.CompAllLineNormals();
}
// Output Cell vertices.
void NewCell_3D ::OutputVertices()
{
	std::cout << "Output cell " << cell_index << " vertices: \n" ;
	std::cout << "vertices.size(): " << vertices.size() << "\n";
	for (unsigned int i = 0; i<vertices.size(); ++i)
		std::cout << "vertices[" << i << "]: "<< vertices[i] << "\n";
}

// Create a vector of vertices with non-repeated vertices in order to calculate the centroid of the cell.
void NewCell_3D ::OrganizeVertices()
{
	for (int face_it = 0; face_it < number_of_faces; ++face_it)
	{
		for (int line_it = 0; line_it < Obj_VectorNewFace[face_it].number_of_lines; ++line_it)
		{
			if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0)
			== vertices.end())
			{
				vertices.push_back(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0);
				vertices_index++;
			}

			if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1)
			== vertices.end())
			{
				vertices.push_back(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1);
				vertices_index++;
			}
		}
	}
	number_of_vertices = vertices.size();
}

// Compute the geometric center of the cell.
void NewCell_3D::CompCellCentroid()
{
	OrganizeVertices();
//	assert(vertices.size()>2);
	double X_centroid = 0;
	double Y_centroid = 0;
	double Z_centroid = 0;
	for (unsigned int i = 0; i<vertices.size(); ++i)
	{
		X_centroid += vertices[i][0];
		Y_centroid += vertices[i][1];
		Z_centroid += vertices[i][2];
	}
	cell_centroid[0] = X_centroid/vertices.size();
	cell_centroid[1] = Y_centroid/vertices.size();
	cell_centroid[2] = Z_centroid/vertices.size();
}
// Reorder all line vertices in counter clockwise order, sort vertices in consecutive arbitrary order and compute face area.
// Need to be called AFTER setting up normal vectors!
void NewCell_3D::ReorderAllVertices()
{
	for (int face_it = 0; face_it < number_of_faces; ++face_it)
	{
		Obj_VectorNewFace[face_it].ReorderAllLineVertices();
		Obj_VectorNewFace[face_it].SortVertices();
		Obj_VectorNewFace[face_it].CompFaceArea();
	}
}
// Checked with volInt from Mirtich and works.
// The first cell of a grid with total refinement = 3 was evaluated by inputting its vertices in
// tetra_test @ Documents/volumeIntegration and its volume was calculated by Mirtich's code (T1)
void NewCell_3D::CompCellVolume()
{
	static const int degree = 1;

	std::vector<double> line_quadrature_weights;
	std::vector<Point<dim-2> > line_quadrature_points;
	unsigned int   line_n_q_points;

	QGauss<dim-2> line_quadrature_formula(2);
	line_n_q_points = line_quadrature_formula.size();
	line_quadrature_weights = line_quadrature_formula.get_weights();
	 line_quadrature_points = line_quadrature_formula.get_points();

	 double H_t_multiplied_sum = 0.0;

	for (int face_it = 0; face_it < number_of_faces; ++face_it)
	{
//		if (Obj_VectorNewFace[face_it].is_boundary_face)
//			std::cout << "IBF \n";

		Point<dim> face_normal = Obj_VectorNewFace[face_it].face_normal_vector;
		double w = Obj_VectorNewFace[face_it].w;
		double gamma_ = Obj_VectorNewFace[face_it].gamma_;
		double alfa = Obj_VectorNewFace[face_it].alfa;
		double beta = Obj_VectorNewFace[face_it].beta;

		for (int line_it = 0; line_it <  Obj_VectorNewFace[face_it].number_of_lines; ++line_it)
		{

			Point<dim-1> line_norm_projection =
					Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].normal_projection;
			double line_length_projection =
					Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].length_projection;
			double real_face_length =
					Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].real_face_length;

			polynomial<double,dim,degree> p0; p0 = 0.0;
			int expon[dim] = {0,0,0};
			p0[p0.term_index(expon)] = 1.;
			assert(p0[0] == 1.);
			std::vector<polynomial<double,dim,degree+1>> p1;
			NPPolynomials::CompF_vector(p0, p1);

			polynomial<double,dim,degree+1> p1_g; p1_g = 0.;

			for (unsigned int coord = 0; coord<dim; ++coord)
				p1_g += p1[coord]*face_normal[coord];

			polynomial<double,dim-1,(degree+1)*2> p1_g_projected =
					NPPolynomials::TransformProjection(p1_g,face_normal/*,line_norm_projection*/,w, alfa, beta, gamma_);

			std::vector<polynomial<double,dim-1,(degree+1)*2+1>> p1_H_vector;
			NPPolynomials::CompF_vector(p1_g_projected, p1_H_vector);
			//////////////////////////////////////////////////////////////////////
//			double H_t_multiplied_sum = 0.0;
			for (unsigned int q_point = 0; q_point<line_n_q_points; ++q_point)
			{
				double alfa_t = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa]+
						(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[alfa] -
								Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[alfa] )*
								line_quadrature_points[q_point](0);

				double beta_t = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta]+
						(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_1_projection[beta] -
								Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X_0_projection[beta] )*
								line_quadrature_points[q_point](0);

//				double alfa_t = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0[alfa]+
//										(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1[alfa] -
//												Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0[alfa] )*
//												line_quadrature_points[q_point](0);
//
//				double beta_t = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0[beta]+
//						(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1[beta] -
//								Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0[beta] )*
//								line_quadrature_points[q_point](0);

				Point<2> point (alfa_t,beta_t);
				Point<2> H_t;
				H_t[X] = NPPolynomials::value(p1_H_vector[X],point);
				H_t[Y] = NPPolynomials::value(p1_H_vector[Y],point);
//				H_t[X] *= line_norm_projection[X]*real_face_length/*line_length_projection*/; // * n_x*1/3
//				H_t[Y] *= line_norm_projection[Y]*real_face_length/*line_length_projection*/; // * n_x*1/3
				H_t[X] *= line_norm_projection[X]*line_length_projection; // * n_x*1/3
				H_t[Y] *= line_norm_projection[Y]*line_length_projection; // * n_x*1/3
				H_t *= (1.0/fabs(face_normal[gamma_]));
				double H_t_multiplied = H_t[X] + H_t[Y];
				H_t_multiplied *= line_quadrature_weights[q_point];
				H_t_multiplied_sum +=H_t_multiplied;
			} // end q points
		} // end loop line
	} // end loop faces
	cell_volume = H_t_multiplied_sum;
}
double NewCell_3D::GetCellVolume()
{
	return cell_volume;
}
